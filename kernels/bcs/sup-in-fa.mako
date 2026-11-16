<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ${c['rho']};
% for i, v in enumerate('uvw'):
    ur[${i + 1}] = (${c['rho']})*(${c[v]});
% endfor
    ur[4] = ${c['Bx']};
    ur[5] = ${c['By']};
    ur[6] = ${c['Bz']};
    ur[7] = ${c['p']}/${c['gamma'] - 1} +
            0.5*(1.0/ur[0])*${pyfr.dot('ur[{i}]', i=(1, 4))} +
            0.5*${pyfr.dot('ur[{i}]', i=(4, 7))};
    ur[8] = ${c['phi']};
</%pyfr:macro>
