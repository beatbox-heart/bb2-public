// Common part of ell_coef.bbs and ell_coefng.bbs: the problem, the solves, the errors.
k_func name=exact when=always pgm={
  def real X (x-1)*hx;
  def real Y (y-1)*hx*gt(ny,0);
  def real Z (z-1)*hx*gt(nz,0);
  def real psi kx*X+ky*Y+kz*Z+ph;
  def real s sin(psi); def real c cos(psi);
  def real k2 kx*kx+ky*ky+kz*kz;
  def real fk fx*kx+fy*ky+fz*kz;
  // isotropic: div(D grad u) = D lap u + grad D . grad u
  def real Di D0*(1+va*sin(bx*X+by*Y));
  def real DiX D0*va*bx*cos(bx*X+by*Y);
  def real DiY D0*va*by*cos(bx*X+by*Y);
  def real fiso -Di*k2*s+(DiX*kx+DiY*ky)*c;
  // anisotropic, constant fibres, Dl=Dl(X): Dt lap u + d/dX(Dl) fx (f.k) c - (Dl-Dt)(f.k)^2 s
  def real Dl Dl0*(1+va*sin(bx*X));
  def real DlX Dl0*va*bx*cos(bx*X);
  def real fani -Dt0*k2*s+DlX*fx*fk*c-(Dl-Dt0)*fk*fk*s;
  u[ua]=s;
  u[f]=iso*fiso+(1-iso)*fani;
  u[p]=u[ua];
  u[g]=u[ua];
};
k_func name=zero when=always x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 pgm={u[p]=0; u[g]=0;};

// Optionally replace the RHS by the diff device's discretization of div(D grad u_exact)
diff when=rhs v0=[ua] v1=[f] [5] [6] hx=hx d0=d0 d1=d1;

// The solves
elliptic when=always x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 v0=[f] v1=[p]
  [5] [6] hx=hx d0=d0 d1=d1 tolerance=tol maxiter=1
  upper_level=ul smoother=GS FMGlevel=0 vcycles=vcyc preiter=2 postiter=2
  debug=[0]_p.dbg;
elliptic when=gs x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 v0=[f] v1=[g]
  [5] [6] hx=hx d0=d0 d1=d1 tolerance=tol maxiter=10000000
  upper_level=-1
  debug=[0]_g.dbg;

// Norms of the errors over the interior
def real Eu; def real Eg; def real Epg;
k_func name=errp x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 pgm={u[d]=max(u[p]-u[ua],u[ua]-u[p])};
reduce x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 v0=[d] v1=[d] operation="max" result=Eu;
k_func name=errg x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 pgm={u[d]=max(u[g]-u[ua],u[ua]-u[g])};
reduce x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 v0=[d] v1=[d] operation="max" result=Eg;
k_func name=errpg x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 pgm={u[d]=max(u[p]-u[g],u[g]-u[p])};
reduce x0=xi0 x1=xi1 y0=yi0 y1=yi1 z0=zi0 z1=zi1 v0=[d] v1=[d] operation="max" result=Epg;

k_func when=always nowhere=1 debug=stdout pgm={Eu=Eu; Eg=Eg; Epg=Epg};
stop when=always;
end;
