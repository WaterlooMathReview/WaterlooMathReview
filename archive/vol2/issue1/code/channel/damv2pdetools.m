function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1.5 1 1]);
set(ax,'PlotBoxAspectRatio',[1 1 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-1 1]);
set(ax,'XTick',[ -1.5,...
 -1.45,...
 -1.3999999999999999,...
 -1.3500000000000001,...
 -1.3,...
 -1.25,...
 -1.2,...
 -1.1499999999999999,...
 -1.1000000000000001,...
 -1.05,...
 -1,...
 -0.94999999999999996,...
 -0.89999999999999991,...
 -0.84999999999999998,...
 -0.79999999999999993,...
 -0.75,...
 -0.69999999999999996,...
 -0.64999999999999991,...
 -0.59999999999999998,...
 -0.54999999999999993,...
 -0.5,...
 -0.44999999999999996,...
 -0.39999999999999991,...
 -0.34999999999999987,...
 -0.29999999999999982,...
 -0.25,...
 -0.19999999999999996,...
 -0.14999999999999991,...
 -0.099999999999999867,...
 -0.049999999999999822,...
 0,...
 0.049999999999999822,...
 0.099999999999999867,...
 0.14999999999999991,...
 0.19999999999999996,...
 0.25,...
 0.29999999999999982,...
 0.34999999999999987,...
 0.39999999999999991,...
 0.44999999999999996,...
 0.5,...
 0.54999999999999993,...
 0.59999999999999998,...
 0.64999999999999991,...
 0.69999999999999996,...
 0.75,...
 0.79999999999999993,...
 0.84999999999999998,...
 0.89999999999999991,...
 0.94999999999999996,...
 1,...
 1.05,...
 1.1000000000000001,...
 1.1499999999999999,...
 1.2,...
 1.25,...
 1.3,...
 1.3500000000000001,...
 1.3999999999999999,...
 1.45,...
 1.5,...
]);
set(ax,'YTick',[ -1,...
 -0.94999999999999996,...
 -0.90000000000000002,...
 -0.84999999999999998,...
 -0.80000000000000004,...
 -0.75,...
 -0.69999999999999996,...
 -0.64999999999999991,...
 -0.59999999999999998,...
 -0.55000000000000004,...
 -0.5,...
 -0.44999999999999996,...
 -0.39999999999999991,...
 -0.34999999999999998,...
 -0.29999999999999993,...
 -0.25,...
 -0.19999999999999996,...
 -0.14999999999999991,...
 -0.099999999999999978,...
 -0.049999999999999933,...
 0,...
 0.049999999999999933,...
 0.099999999999999978,...
 0.14999999999999991,...
 0.19999999999999996,...
 0.25,...
 0.29999999999999993,...
 0.34999999999999998,...
 0.39999999999999991,...
 0.44999999999999996,...
 0.5,...
 0.55000000000000004,...
 0.59999999999999998,...
 0.64999999999999991,...
 0.69999999999999996,...
 0.75,...
 0.80000000000000004,...
 0.84999999999999998,...
 0.90000000000000002,...
 0.94999999999999996,...
 1,...
]);
pdetool('gridon','on');

% Geometry description:
pdepoly([ 0.5,...
 -0.24999999999999978,...
 -0.19999999999999996,...
 -0.099999999999999867,...
 2.2204460492503131e-016,...
 0.14999999999999991,...
 0.25,...
 0.34999999999999987,...
 0.5,...
],...
[ -0.69999999999999996,...
 -0.69999999999999996,...
 -0.59999999999999998,...
 -0.44999999999999996,...
 -0.34999999999999998,...
 -0.25,...
 -0.20000000000000007,...
 -0.15000000000000002,...
 -0.10000000000000009,...
],...
 'P1');
pdepoly([ 0.49999999999999978,...
 -0.25,...
 -0.20000000000000018,...
 -0.10000000000000009,...
 0,...
 0.14999999999999969,...
 0.24999999999999978,...
 0.34999999999999964,...
 0.49999999999999978,...
],...
[ 0.69999999999999951,...
 0.69999999999999951,...
 0.59999999999999953,...
 0.44999999999999962,...
 0.34999999999999953,...
 0.24999999999999956,...
 0.19999999999999962,...
 0.14999999999999958,...
 0.099999999999999534,...
],...
 'P2');
pderect([1.2 -1.45 0.69999999999999996 -0.69999999999999996],'R1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-P2-P1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(22,...
'dir',...
1,...
'1',...
'0')
pdesetbd(21,...
'dir',...
1,...
'1',...
'0')
pdesetbd(20,...
'dir',...
1,...
'1',...
'0')
pdesetbd(19,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(17,...
'dir',...
1,...
'1',...
'0')
pdesetbd(16,...
'dir',...
1,...
'1',...
'0')
pdesetbd(15,...
'dir',...
1,...
'1',...
'0')
pdesetbd(14,...
'dir',...
1,...
'1',...
'0')
pdesetbd(13,...
'dir',...
1,...
'1',...
'0')
pdesetbd(12,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'0')
pdesetbd(8,...
'dir',...
1,...
'1',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(6,...
'dir',...
1,...
'1',...
'0')
pdesetbd(5,...
'dir',...
1,...
'1',...
'0')
pdesetbd(4,...
'dir',...
1,...
'1',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'0')
pdesetbd(2,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
str2mat('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');
