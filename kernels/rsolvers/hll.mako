<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.mhd.kernels.flux'/>

## HLL solver for GLM-MHD. The wave speed bounds use the fast magnetosonic
## speed and are widened to include the divergence cleaning waves at +/- vref,
## so sl < 0 < sr always holds and no upwind shortcut is needed.
<%pyfr:macro name='rsolve' params='ul, ur, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[3], vr[3];
    fpdtype_t pl, pr;
    fpdtype_t nf_fl, nf_fr;

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

    // Estimate the left and right wave speeds
    fpdtype_t sl = fmin(fmin(nvl - cfl, nvr - cfr), -${c['vref']});
    fpdtype_t sr = fmax(fmax(nvl + cfl, nvr + cfr), ${c['vref']});
    fpdtype_t rcpsrsl = 1/(sr - sl);

    // Output
% for i in range(nvars):
    nf_fl = ${' + '.join(f'n[{j}]*fl[{j}][{i}]' for j in range(ndims))};
    nf_fr = ${' + '.join(f'n[{j}]*fr[{j}][{i}]' for j in range(ndims))};
    nf[${i}] = (sr*nf_fl - sl*nf_fr + sl*sr*(ur[${i}] - ul[${i}]))*rcpsrsl;
% endfor
</%pyfr:macro>
