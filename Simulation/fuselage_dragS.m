function [Cdf, S] = fuselage_dragS()
    % Physical quanities of the helicopter
    fus_width = 0.99; %m
    fus_length  = 13.77; %m
    fus_height = 3.11; %m
    wings_width = 3.56-0.99; %m the wing width is the total width of the wings minus the fuselage width
    wings_thickness = 0.1; %m an approximation of the wing thickness
    
    Ld = fus_length/fus_width;
    
    FF = 1 + 2.2/(Ld*1.5)^1.5 + 3.8/(Ld)^3;  %form factor
    
    Re = 1e6; %Reynolds number 
    
    Cf = (1/(3.46*log10(Re)-5.6))^2; %skin friction coefficient
    
    S = (fus_length*fus_width + fus_height*fus_length + fus_height*fus_width) *2* 2/3; %2/3 correction for the tail getting narrower
    
    Fus_area = fus_width*fus_height + wings_width*wings_thickness;
    
    Cdf = Cf *FF *S/Fus_area;
end