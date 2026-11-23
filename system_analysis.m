% A matlab program to plot the theoretical response of the building in the
% 1A vibration lab.
% Based on code written by Penny Cox, now maintained by Aidan Reilly.
%  Tidied up by Jim Woodhouse October 2012
% Code applied by Yitian Chen November 2026


close all
clear all

%  USER MUST SET TWO VALUES TO DETERMINE PLOT OPTIONS

% Set value for plottype as follows:
%       1 -> use 'ezplot' to plot a graph on normal axis (in rad/s)
%       2 -> produce a semilog plot (in Hz)
%       3 -> produce a graph on linear axis using 'plot' (in Hz)
plottype=1;

% Set value to be positive to produce plots of the modeshapes
modeshape_visualisation = 0;

m = 1.83; % mass of one floor
L = 0.2; % length
N = 3; % number of degrees of freedom
b = 0.08; % width
E = 210E9; % Young's Modulus
d = 0.001; % thickness
I = b*d*d*d/12; % second moment of area
k = (24*E*I)/(L*L*L); % static stiffness for each floor

M = m*eye(N); % create the mass matrix
K = k*[2 -1 0;-1 2 -1;0 -1 1]; % create the stiffness matrix

% Modified the M and K matrices above to include the effect of mass
% absorbers on each floor

% Natural frequencies of frame without absorbers sorted by lowest to
% highest
w0 = sort(sqrt(eig(K, M)));

% Absorber mass at each floor
m_abs = [0.1, 0.2, 0.1];

% Absorber property matrix of rows, following the convention [floor, absorber mass,
% stiffness]. 

absorbers = [ 3, m_abs(3), m_abs(3) * w0(1)^2; % floor 3 tuned to mode 1
    2, m_abs(2), m_abs(2) * w0(3)^2; % floor 2 tuned to mode 3
    1, m_abs(1), m_abs(1) * w0(2)^2;]; % floor 1 tuned to mode 2

% Add the absorber properties to M and K matrices one at a time

for j = 1:size(absorbers, 1)

    floor = absorbers(j, 1);
    absorber_mass = absorbers(j, 2);
    absorber_stiffness = absorbers(j, 3);

    i = size(M, 1) + 1; % index of absorber DOF

    M(i, i) = absorber_mass; 
    K(i, i) = absorber_stiffness;
    K(floor, floor) = K(floor, floor) + absorber_stiffness;
    K(i, floor) = -absorber_stiffness;
    K(floor, i) = -absorber_stiffness;

end

N = size(M, 1); % should be 6 DOF with 3 floors and 3 absorbers

F = zeros(N, 1); % force vector

F(1) = 1; % unit force applied to floor 1

[V,D] = eig(K,M);
syms w;

for imode=1:N
  freqs(imode) = sqrt(D(imode,imode));
end

%  Print natural frequencies and mode vectors in command window
hertz = freqs/(2*pi)
modeshapes = V

B = K - ((w*w)*M); 
% harmonic solution for unit force at floor 1
disp = (inv(B))*F;

%start of ezplot section
if (plottype == 1)
  hold on
  
  ifloor=1;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line'),'Color','k')
  ifloor=2;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','g')
  ifloor=3;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','r')
  
  set(findobj('Type','line','Color','k'),'Color','b')
  
  set(findobj('Type','line'),'LineStyle','-')
end

% Calculate frequency response functions
all_disp = [];
for w = 1:130;
  B = K - ((w*w)*M); 
  % harmonic solution for unit force at floor 1
  disp = (inv(B))*F;
  all_disp = [all_disp disp];
end
w = 1:130;

% Log plot
if (plottype == 2)
  semilogy((w./(2*pi)),abs(all_disp),'-');

% Linear plot
elseif (plottype == 3)
  plot((w./(2*pi)),(all_disp),'-');
end

% Plot modeshapes

if (modeshape_visualisation > 0 )
  V = [0 0 0; V];
  V_ = V + 0.25;
  V = V - 0.25;
  for imode=1:3
    figure
    axis([-5 5 0 3.5])
    title1 = ['Mode ' int2str(imode)];
    title(title1)
    hold on
    plot((V(:,imode)),([0 1 2 3]))
    plot([0 0 0 0],[0 1 2 3],'k')
    plot((V_(:,imode)),([0 1 2 3]))
    for jmode=1:3
      plot([V(jmode+1,imode) V_(jmode+1,imode)],[jmode jmode])
    end
  end
end