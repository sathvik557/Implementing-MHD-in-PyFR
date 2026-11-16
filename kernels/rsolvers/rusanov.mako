<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mhd.kernels.flux'/>

## Rusanov (local Lax-Friedrichs) solver for GLM-MHD.
## Maximum wave speed: the larger of |v_n| + c_fast (left or right state) and
## the divergence cleaning speed vref.
<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[3], vr[3];
    fpdtype_t pl, pr;

    ${pyfr.expand('inviscid_flux', 'ul', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'fr', 'pr', 'vr')};

    // Normal velocities and normal magnetic fields
    fpdtype_t nvl = ${pyfr.dot('n[{i}]', 'vl[{i}]', i=ndims)};
    fpdtype_t nvr = ${pyfr.dot('n[{i}]', 'vr[{i}]', i=ndims)};
    fpdtype_t nbl = ${pyfr.dot('n[{i}]', 'ul[{i} + 4]', i=ndims)};
    fpdtype_t nbr = ${pyfr.dot('n[{i}]', 'ur[{i} + 4]', i=ndims)};

    // Fast magnetosonic speeds
    fpdtype_t a2l = ${c['gamma']}*pl/ul[0];
    fpdtype_t a2r = ${c['gamma']}*pr/ur[0];
    fpdtype_t b2l = ${pyfr.dot('ul[{i}]', i=(4, 7))}/ul[0];
    fpdtype_t b2r = ${pyfr.dot('ur[{i}]', i=(4, 7))}/ur[0];
    fpdtype_t cfl = sqrt(0.5*(a2l + b2l + sqrt(fmax((a2l + b2l)*(a2l + b2l)
                                                  - 4*a2l*nbl*nbl/ul[0], 0))));
    fpdtype_t cfr = sqrt(0.5*(a2r + b2r + sqrt(fmax((a2r + b2r)*(a2r + b2r)
                                                  - 4*a2r*nbr*nbr/ur[0], 0))));

    // Estimate the maximum wave speed
    fpdtype_t lam = fmax(fmax(fabs(nvl) + cfl, fabs(nvr) + cfr),
                         ${c['vref']});

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join(f'n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 for j in range(ndims))})
             + 0.5*lam*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
