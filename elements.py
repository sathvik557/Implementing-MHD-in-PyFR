from pyfr.solvers.baseadvec import BaseAdvectionElements


class BaseMHDElements:
    # 2.5D: always three velocity and three B components, even on a 2D mesh
    @staticmethod
    def privars(ndims, cfg):
        return ['rho', 'u', 'v', 'w', 'Bx', 'By', 'Bz', 'p']

    @staticmethod
    def convars(ndims, cfg):
        return ['rho', 'rhou', 'rhov', 'rhow', 'Bx', 'By', 'Bz', 'E']

    dualcoeffs = convars

    @staticmethod
    def visvars(ndims, cfg):
        return {
            'density': ['rho'],
            'velocity': ['u', 'v', 'w'],
            'magnetic-field': ['Bx', 'By', 'Bz'],
            'pressure': ['p']
        }

    @staticmethod
    def pri_to_con(pris, cfg):
        rho, p = pris[0], pris[-1]
        vs, Bs = pris[1:4], pris[4:7]

        gamma = cfg.getfloat('constants', 'gamma')
        mu0 = cfg.getfloat('constants', 'mu0', 1.0)

        rhovs = [rho*v for v in vs]
        E = (p/(gamma - 1) + 0.5*rho*sum(v*v for v in vs)
             + sum(b*b for b in Bs)/(2*mu0))

        return [rho, *rhovs, *Bs, E]

    @staticmethod
    def con_to_pri(cons, cfg):
        rho, E = cons[0], cons[-1]
        rhovs, Bs = cons[1:4], cons[4:7]

        gamma = cfg.getfloat('constants', 'gamma')
        mu0 = cfg.getfloat('constants', 'mu0', 1.0)

        vs = [rhov/rho for rhov in rhovs]
        p = (gamma - 1)*(E - 0.5*rho*sum(v*v for v in vs)
                         - sum(b*b for b in Bs)/(2*mu0))

        return [rho, *vs, *Bs, p]

    @staticmethod
    def diff_con_to_pri(cons, diff_cons, cfg):
        rho = cons[0]
        rhovs, Bs = cons[1:4], cons[4:7]
        diff_rho, diff_E = diff_cons[0], diff_cons[-1]
        diff_rhovs, diff_Bs = diff_cons[1:4], diff_cons[4:7]

        gamma = cfg.getfloat('constants', 'gamma')
        mu0 = cfg.getfloat('constants', 'mu0', 1.0)

        vs = [rhov/rho for rhov in rhovs]
        diff_vs = [(drv - v*diff_rho)/rho for drv, v in zip(diff_rhovs, vs)]

        # Kinetic part is the same as Euler; magnetic part is B.dB/mu0
        diff_p = diff_E - 0.5*(sum(v*drv for v, drv in zip(vs, diff_rhovs)) +
                               sum(rv*dv for rv, dv in zip(rhovs, diff_vs)))
        diff_p -= sum(b*db for b, db in zip(Bs, diff_Bs))/mu0
        diff_p *= gamma - 1

        return [diff_rho, *diff_vs, *diff_Bs, diff_p]

    @staticmethod
    def validate_formulation(ctrl):
        shock_capturing = ctrl.cfg.get('solver', 'shock-capturing', 'none')
        if shock_capturing == 'entropy-filter':
            raise ValueError('Entropy filtering is not implemented for MHD.')


class MHDElements(BaseMHDElements, BaseAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux kernels
        self._be.pointwise.register('pyfr.solvers.mhd.kernels.tflux')

        # Template parameters for the flux kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs
        }

        # Helpers
        tdisf = []
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat
        slicedk = self._make_sliced_kernel

        if c in r and 'flux' not in self.antialias:
            tdisf.append(lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'curved'},
                dims=[self.nupts, r[c]], u=s(self.scal_upts[uin], c),
                f=s(self._vect_upts, c), smats=self.curved_smat_at('upts')
            ))
        elif c in r:
            tdisf.append(lambda: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'curved'},
                dims=[self.nqpts, r[c]], u=s(self._scal_qpts, c),
                f=s(self._vect_qpts, c), smats=self.curved_smat_at('qpts')
            ))

        if l in r and 'flux' not in self.antialias:
            tdisf.append(lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'linear'},
                dims=[self.nupts, r[l]], u=s(self.scal_upts[uin], l),
                f=s(self._vect_upts, l), verts=self.ploc_at('linspts', l),
                upts=self.upts
            ))
        elif l in r:
            tdisf.append(lambda: self._be.kernel(
                'tflux', tplargs=tplargs | {'ktype': 'linear'},
                dims=[self.nqpts, r[l]], u=s(self._scal_qpts, l),
                f=s(self._vect_qpts, l), verts=self.ploc_at('linspts', l),
                upts=self.qpts
            ))

        if 'flux' not in self.antialias:
            self.kernels['tdisf'] = lambda uin: slicedk(k(uin) for k in tdisf)
        else:
            self.kernels['tdisf'] = lambda: slicedk(k() for k in tdisf)
