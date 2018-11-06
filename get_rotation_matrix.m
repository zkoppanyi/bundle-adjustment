function R = get_rotation_matrix(omega, phi, kappa)

    if nargin == 1
        x = omega;
        omega = x(1);
        phi = x(2);
        kappa = x(3);
    end
    
    R = zeros(3,3);
    R(1,1) = cos(phi)*cos(kappa);
    R(1,2) = -cos(phi)*sin(kappa);
    R(1,3) = sin(phi);
    R(2,1) = cos(omega)*sin(kappa)+sin(omega)*sin(phi)*cos(kappa);
    R(2,2) = cos(omega)*cos(kappa)-sin(omega)*sin(phi)*sin(kappa);
    R(2,3) = -sin(omega)*cos(phi);
    R(3,1) = sin(omega)*sin(kappa)-cos(omega)*sin(phi)*cos(kappa);
    R(3,2) = sin(omega)*cos(kappa)+cos(omega)*sin(phi)*sin(kappa);
    R(3,3)  = cos(omega)*cos(phi);

%     Rx = [1 0 0; 0 cos(omega) -sin(omega); 0 sin(omega) cos(omega)];
%     Ry = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];
%     Rz = [cos(kappa) -sin(kappa) 0; sin(kappa) cos(kappa) 0; 0 0 1];
%     R = Rz * Ry * Rx;

