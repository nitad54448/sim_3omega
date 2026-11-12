// Uniform heat flux over a strip of width 2b; temperature averaged over that strip.
// Convention: Re[T exp(i 2 omega t)], current specified as a peak amplitude.
const EULER_GAMMA = 0.5772156649015329;
const CM = {
 add:(a,b)=>({re:a.re+b.re,im:a.im+b.im}),
 sub:(a,b)=>({re:a.re-b.re,im:a.im-b.im}),
 scale:(a,s)=>({re:a.re*s,im:a.im*s}),
 mul:(a,b)=>({re:a.re*b.re-a.im*b.im,im:a.re*b.im+a.im*b.re}),
 div:(a,b)=>{const t=Math.max(Math.abs(b.re),Math.abs(b.im));if(!(t>0)||!Number.isFinite(t))throw Error('Invalid thermal denominator.');const r=b.re/t,i=b.im/t,d=r*r+i*i;return {re:(a.re/t*r+a.im/t*i)/d,im:(a.im/t*r-a.re/t*i)/d};},
 abs:a=>Math.hypot(a.re,a.im),
 sqrtPositive:(x,y)=>{const r=Math.hypot(x,y),re=Math.sqrt(r/2+x/2);return {re,im:y/(2*re)};},
 tanh:z=>{
  // exp(-2z) avoids overflowing hyperbolic functions for thick substrates.
  if(z.re>20)return {re:1,im:0};
  if(CM.abs(z)<1e-4){const z2=CM.mul(z,z);return CM.mul(z,CM.add({re:1,im:0},CM.add(CM.scale(z2,-1/3),CM.scale(CM.mul(z2,z2),2/15))));}
  const e=Math.exp(-2*z.re),c=Math.cos(2*z.im),s=Math.sin(2*z.im);
  return CM.div({re:-Math.expm1(-2*z.re)+e*(1-c),im:e*s},{re:1+e*c,im:-e*s});
 }
};
function adaptiveComplexSimpson(fn,a,b,tolerance=1e-10) {
 const evalAt=x=>{const v=fn(x);if(!Number.isFinite(v.re)||!Number.isFinite(v.im))throw Error('Nonfinite thermal integral. Check parameter magnitudes.');return v;};
 const simpson=(a,b,fa,fm,fb)=>CM.scale(CM.add(CM.add(fa,CM.scale(fm,4)),fb),(b-a)/6);
 const fa=evalAt(a),fb=evalAt(b),fm=evalAt((a+b)/2),whole=simpson(a,b,fa,fm,fb);
 const refine=(a,b,fa,fm,fb,whole,tol,depth)=>{
  const m=(a+b)/2,l=evalAt((a+m)/2),r=evalAt((m+b)/2);
  const left=simpson(a,m,fa,l,fm),right=simpson(m,b,fm,r,fb),sum=CM.add(left,right),delta=CM.sub(sum,whole),error=CM.abs(delta)/15;
  if(error<=tol)return {value:CM.add(sum,CM.scale(delta,1/15)),error};
  if(depth===0)throw Error('Integration did not converge. Parameters are outside the supported numerical range.');
  const A=refine(a,m,fa,l,fm,left,tol/2,depth-1),B=refine(m,b,fm,r,fb,right,tol/2,depth-1);
  return {value:CM.add(A.value,B.value),error:A.error+B.error};
 };
 return refine(a,b,fa,fm,fb,whole,tolerance,26);
}
function thermalIntegral(f,b,D,thickness,boundary='infinite') {
 if(![f,b,D,thickness].every(v=>Number.isFinite(v)&&v>0))throw Error('Frequency, dimensions and diffusivity must be finite and positive.');
 if(!['infinite','isothermal'].includes(boundary))throw Error('Unknown substrate boundary condition.');
 const a=4*Math.PI*f*b*b/D,r=thickness/b;
 if(!Number.isFinite(a)||a<=0||!Number.isFinite(r))throw Error('Parameter magnitudes exceed the supported numerical range.');
 const kernel=x=>{
  const g=CM.sqrtPositive(x*x,a),sinc=x===0?1:Math.sin(x)/x;
  let impedance;
  if(boundary==='isothermal'){
   const z=CM.scale(g,r);
   // tanh(g r)/g -> r as g r -> 0; retain the complex series.
   if(CM.abs(z)<1e-4){const z2=CM.mul(z,z);impedance=CM.scale(CM.add({re:1,im:0},CM.add(CM.scale(z2,-1/3),CM.scale(CM.mul(z2,z2),2/15))),r);}
   else impedance=CM.div(CM.tanh(z),g);
  }else impedance=CM.div({re:1,im:0},g);
  return CM.scale(impedance,sinc*sinc);
 };
 let value={re:0,im:0},error=0;
 const addInterval=(lo,hi)=>{const v=adaptiveComplexSimpson(kernel,lo,hi);value=CM.add(value,v.value);error+=v.error;};
 // Resolve the low-x scale sqrt(a), then sample every sinc lobe separately.
 let x=Math.min(1,Math.sqrt(a));
 addInterval(0,x);
 while(x<1){const next=Math.min(1,2*x);if(next===x)throw Error('Unresolvable thermal length scale.');addInterval(x,next);x=next;}
 let cutoff=256;
 for(;;){
  while(x<cutoff){const next=Math.min(cutoff,x+Math.PI/2);addInterval(x,next);x=next;}
  // For x>=X: sinc²(x)<=1/x², |g|>=x and |tanh(g r)|<=coth(X r).
  const boundFactor=boundary==='isothermal'?1/Math.tanh(cutoff*r):1;
  const tailBound=boundFactor/(2*cutoff*cutoff),estimatedError=error+tailBound;
  if(estimatedError<=2e-5*CM.abs(value))return {...value,estimatedError,cutoff};
  if(cutoff>=8192)throw Error('Integral tail tolerance not reached; parameter combination is outside the supported numerical range.');
  cutoff*=2;
 }
}
