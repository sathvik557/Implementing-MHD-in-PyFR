from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)


class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        self.p_min = self.cfg.getfloat('solver-interfaces', 'p-min',
                                       5*self._be.fpdtype_eps)

        # Divergence cleaning speed (V_ref,HDC in Dhib et al. 2026). It has
        # no safe default: it should be about the fastest wave in the problem
        if 'vref' not in self.c:
            raise ValueError('MHD needs "vref" (the divergence cleaning '
                             'speed) in the [constants] section')

        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver, c=self.c, p_min=self.p_min)


class MHDIntInters(TplargsMixin, BaseAdvectionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mhd.kernels.intcflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class MHDMPIInters(TplargsMixin, BaseAdvectionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mhd.kernels.mpicflux')

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', self._tplargs, dims=[self.ninterfpts],
            ul=self._scal_lhs, ur=self._scal_rhs, nl=self._pnorm_lhs
        )


class MHDBaseBCInters(TplargsMixin, BaseAdvectionBCInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.mhd.kernels.bccflux')

        self._tplargs |= dict(bctype=self.type, ninters=self.ninters)

        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs, nl=self._pnorm_lhs,
            **self._external_vals
        )


class MHDSupInflowBCInters(MHDBaseBCInters):
    type = 'sup-in-fa'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w', 'Bx', 'By', 'Bz', 'phi'], lhs
        )


class MHDSupOutflowBCInters(MHDBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class MHDSlpAdiaWallBCInters(MHDBaseBCInters):
    type = 'slp-adia-wall'
