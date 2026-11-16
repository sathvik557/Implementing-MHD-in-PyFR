<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

## Slip wall: mirror the normal momentum only. The magnetic field, energy and
## cleaning variable are copied from the inside state; whether that is the
## right treatment for a perfectly conducting wall is a modelling choice.
<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    fpdtype_t nor = ${' + '.join(f'ul[{i + 1}]*nl[{i}]' for i in range(ndims))};
    ur[0] = ul[0];
% for i in range(3):
% if i < ndims:
    ur[${i + 1}] = ul[${i + 1}] - 2*nor*nl[${i}];
% else:
    ur[${i + 1}] = ul[${i + 1}];
% endif
% endfor
% for i in range(4, nvars):
    ur[${i}] = ul[${i}];
% endfor
</%pyfr:macro>
