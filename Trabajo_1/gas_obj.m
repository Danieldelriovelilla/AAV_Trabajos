classdef gas_obj < handle

   properties (Access = public)
       
       Mm       % Masa molecular [g/mol]
       gamma    % relacion calores
       R
       cp
       cv
       a
       
       T
       P
       rho
       e
       h
       
       T0
       P0
       rho0
       h0
       
   end
   
   methods %(Static)
        
        % Constructor
        function [] = gas_TP(obj,Mm, gamma, P, T)    
            obj.Mm = Mm;
            obj.gamma = gamma;
            
            obj.R = 8314/Mm;
            obj.cp = gamma*obj.R/(gamma-1);
            obj.cv = obj.R/(gamma-1);
            obj.a = sqrt(obj.gamma*obj.R*obj.T);
            
            obj.P = P;
            obj.T = T;
            obj.rho = P/(obj.R*T);
            
            obj.e = obj.cv*T;
            obj.h = obj.cp*T;
            
            obj.a = sqrt(gamma*obj.R*T);
        end
        
        function [T0, P0, rho0, h0] = remanso(obj, M, T, P, rho)
            
            gamma = obj.gamma;
            cp = obj.cp;
            
            T0 = T*( 1 + (gamma-1)/2*M^2 );
            P0 = P*( ( 1 + (gamma-1)/2*M^2 )^(gamma/(gamma-1)) );
            rho0 = rho*( ( 1 + (gamma-1)/2*M^2 )^(1/(gamma-1)) );
            h0 = T0*cp;
            %P0 = P*(T0/obj.T)^((gamma-1)/gamma);
            %rho0 = rho*(T0/T)^(gamma-1);
            

        end
        
        function [M2, T2, P2, rho2, h2] = NOCH(obj, M1, T1, P1, rho1, h1)
        
            gamma = obj.gamma;
            
            M2 = sqrt( ( 1 + ((gamma-1)/2)*M1^2 )/...
                ( gamma*M1^2 - (gamma-1)/2 ) );
            
            T2 = T1*( 1 + 2*gamma/(gamma+1)*(M1^2-1) )*...
                ( (2+(gamma-1)*M1^2)/((gamma+1)*M1^2) );
            
            P2 = P1*( 1 + 2*gamma/(gamma+1)*(M1^2-1) );
            
            rho2 = rho1*( ( (gamma+1)*M1^2 )/...
                ( 2 + (gamma-1)*M1^2 ) );
            
            
            h2 = h1*(T2/T1);
            
        end
        
   end
end