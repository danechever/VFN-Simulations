%% DM parameters 

dm.Nact = 34;	% # of actuators across DM array
dm.dm_spacing = 400e-6;	% User defined actuator pitch
dm.xc = (dm.Nact/2 - 1/2);	% x-center location of DM surface [actuator widths]
dm.yc = (dm.Nact/2 - 1/2);	% y-center location of DM surface [actuator widths]
dm.xtilt = 0; % for foreshortening. angle of rotation about x-axis [degrees]
dm.ytilt = 0; % for foreshortening. angle of rotation about y-axis [degrees]
dm.zrot = 0;  % clocking of DM surface [degrees]
dm.Vmin = 0;  % Min voltage (V)
dm.Vmax = 220;% Max voltage (V)
dm.VtoH = 12.4e-9*ones(dm.Nact);  % gains of all actuators [nm/V of free stroke]
dm.biasMap = 110*ones(dm.Nact); %zeros(dm.Nact); % Bias map (V)
dm.pinned = []; % Pinned actuators 
dm.Vpinned = [];% Voltage of pinned acts 
dm.tied = []; % Tied actuators 
dm.flagNbrRule = false; % Whether to use a neighbor rule 
dm.fitType = 'linear'; % Displacement versus voltage model 

p1 = 0.035e-9;
p2 = 4.58e-9;
p3 = -49.9e-9;

% Influence function 
dm.inf_sign = '+';
dm.inf_fn = 'influence_BMC_2kDM_400micron_res10.fits';

dm.centering = 'pixel'; % Use pixel unless you have a reason not to 

