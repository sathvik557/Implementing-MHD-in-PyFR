<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## Ideal MHD flux with GLM divergence cleaning (Dhib et al. 2026, Eq. 22).
##
## State layout (see convars in elements.py):
##   s = [rho, rhou, rhov, rhow, Bx, By, Bz, E, phi]     (nvars = 9)
## The velocity and magnetic field always have three components (2.5D), but
## there is one flux direction per spatial dimension (ndims).
<%pyfr:macro name='inviscid_flux' params='s, f, p, v'>
    fpdtype_t invrho = 1.0/s[0], E = s[7];

    // Momenta, velocities and magnetic field (always three components)
    fpdtype_t rhov[3], B[3];
% for i in range(3):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
    B[${i}] = s[${i + 4}];
% endfor

    // Magnetic pressure and v.B
    fpdtype_t pm = 0.5*${pyfr.dot('B[{i}]', i=3)};
    fpdtype_t vB = ${pyfr.dot('v[{i}]', 'B[{i}]', i=3)};

    // Gas pressure
    p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=3)} - pm);

    // Flux in each spatial direction
% for i in range(ndims):
    // Direction ${i}: density, energy and cleaning-variable fluxes
    f[${i}][0] = rhov[${i}];
    f[${i}][7] = (E + p + pm)*v[${i}] - B[${i}]*vB;
    f[${i}][8] = ${c['vref']**2}*B[${i}];
% endfor

    // Momentum fluxes: rho v_i v_j + (p + |B|^2/2) delta_ij - B_i B_j
% for i, j in pyfr.ndrange(ndims, 3):
    f[${i}][${j + 1}] = rhov[${i}]*v[${j}] - B[${i}]*B[${j}]${' + p + pm' if i == j else ''};
% endfor

    // Induction fluxes: v_i B_j - B_i v_j + phi delta_ij
% for i, j in pyfr.ndrange(ndims, 3):
    f[${i}][${j + 4}] = v[${i}]*B[${j}] - B[${i}]*v[${j}]${' + s[8]' if i == j else ''};
% endfor
</%pyfr:macro>
