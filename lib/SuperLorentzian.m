%% Super Lorentzian, interpolated as in Gloor et al, 2008
% return LUT interpolated from -1kHz to +1KHz
% Units of G are us (microseconds)

function [ff,G] = SuperLorentzian(T2b)

%%% define frequency range
n=512;
ff = linspace(-30e3,30e3,n);

%%% compute G for this range
G = zeros([n 1]);

for ii=1:n
    G(ii) = SL(ff(ii));
end

%%% interpolate
po = find(abs(ff)<1.5e3); % points to interpolate
pu = find((abs(ff)>1.5e3)&(abs(ff)<2e3)); % points to use

Gi = spline(ff(pu),G(pu),ff(po));
G(po) = Gi;

G = G*1e6; % us

    function gg = SL(f)
        u = linspace(0,1,500);
        du = u(2)-u(1);
        g = sqrt(2/pi).*T2b./(abs(3*u.^2-1));
        g = g .* exp(-2*(2*pi*f*T2b./(3*u.^2-1)).^2);
        gg = du*sum(g);
    end

end


