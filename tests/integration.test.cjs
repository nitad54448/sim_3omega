const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const root=path.join(__dirname,'..');
const kernel=fs.readFileSync(path.join(root,'thermal-integral.js'),'utf8');
const context=vm.createContext({});vm.runInContext(kernel,context);

// Independent fixed-grid reference after x=sqrt(a)*sinh(u). The Jacobian
// resolves the low-frequency length scale without using the production
// interval subdivision, complex helpers or adaptive error estimator.
function reference(a,n=131072,cutoff=512){
  const scale=Math.sqrt(a),h=Math.asinh(cutoff/scale)/n;
  let re=0,im=0;
  for(let i=0;i<=n;i++){
    const u=i*h,x=scale*Math.sinh(u),x2=x*x,r=Math.hypot(x2,a);
    const denomRe=Math.sqrt((r+x2)/2),denomIm=a/(2*denomRe);
    const sinc=x===0?1:Math.sin(x)/x;
    const weight=(i===0||i===n?1:i%2?4:2)*sinc*sinc*scale*Math.cosh(u)/r;
    re+=weight*denomRe;im-=weight*denomIm;
  }
  return {re:re*h/3,im:im*h/3,tailBound:1/(2*cutoff*cutoff)};
}
test('low-frequency regression resolves the reviewed 0.02 Hz response',()=>{
  const J=context.thermalIntegral(.02,5e-6,1e-6,.0005);
  assert(Math.abs(J.re-6.912)<.001);assert(Math.abs(J.im+.785398)<.00001);
  assert(J.estimatedError/Math.hypot(J.re,J.im)<=2e-5);
});
test('adaptive response and error estimate agree with independently refined transformed quadrature',()=>{
  for(const a of [4*Math.PI*.02*(5e-6)**2/1e-6,.001,1,1000]){
    const coarse=reference(a,65536),fine=reference(a);
    const referenceError=Math.hypot(fine.re-coarse.re,fine.im-coarse.im);
    assert(referenceError<1e-8,`reference convergence at a=${a}`);
    const J=context.thermalIntegral(a/(4*Math.PI),1,1,1);
    assert(Math.hypot(J.re-fine.re,J.im-fine.im)<=J.estimatedError+fine.tailBound+referenceError,`error budget at a=${a}`);
    assert(J.re>0&&J.im<0);
  }
});
test('the low-frequency logarithmic asymptote has slope -1/2 and imaginary part -pi/4',()=>{
  const b=5e-6,D=1e-6,d=.0005,f=.00002;
  const J=context.thermalIntegral(f,b,D,d),K=context.thermalIntegral(f*2,b,D,d);
  const asymptote=-.5*Math.log(4*Math.PI*f*b*b/D)+1.5-.5772156649015329;
  assert(Math.abs(J.re-asymptote)<2e-5);assert(Math.abs(J.im+Math.PI/4)<1e-6);
  assert(Math.abs((K.re-J.re)/Math.log(2)+.5)<1e-6);
});
test('finite isothermal boundary still tends to the infinite response for thick substrates',()=>{
  const infinite=context.thermalIntegral(1,5e-6,1e-6,1,'infinite');
  const finite=context.thermalIntegral(1,5e-6,1e-6,1,'isothermal');
  assert(Math.hypot(finite.re-infinite.re,finite.im-infinite.im)<1e-9);
});
test('invalid inputs and non-finite integrands are refused',()=>{
  for(const f of [0,-1,NaN,Infinity])assert.throws(()=>context.thermalIntegral(f,5e-6,1e-6,.0005),/finite and positive/);
  assert.throws(()=>context.adaptiveComplexSimpson(()=>({re:Infinity,im:0}),0,1),/Nonfinite thermal integral/);
});
function pageHarness(file){
  const html=fs.readFileSync(path.join(root,file),'utf8'),nodes=new Map(),blobs=[];
  for(const m of html.matchAll(/<input\b[^>]*id="([^"]+)"[^>]*value="([^"]*)"/g))nodes.set(m[1],{value:m[2]});
  const node=id=>{if(!nodes.has(id))nodes.set(id,{});const n=nodes.get(id);n.style??={};n.addEventListener??=()=>{};return n;};
  node('boundary-condition').value='infinite';
  node('boundary-condition').options=[{text:'Semi-infinite substrate'}];node('boundary-condition').selectedIndex=0;
  const ctx=vm.createContext({window:{jspdf:{},addEventListener(){}},d3:{select:()=>({})},
    document:{getElementById:node,querySelector:selector=>({textContent:selector}),
      createElement:()=>({style:{},click(){}}),body:{appendChild(){},removeChild(){}}},
    Blob,URL:{createObjectURL:blob=>{blobs.push(blob);return 'blob:test';},revokeObjectURL(){}},
    setTimeout:fn=>{fn();return 1;},console});
  const scripts=[...html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script>/gi)];
  for(const m of scripts){if(m[1].includes('./thermal-integral.js'))vm.runInContext(kernel,ctx);else if(m[2].trim())vm.runInContext(m[2],ctx);}
  ctx.drawPlot=()=>{};
  return {ctx,node,blobs,data:()=>vm.runInContext('fullPlotData',ctx)};
}
for(const file of ['index.html','infinite_substrate.html']){
  test(`${file} uses the shared integrator throughout the sweep and clears failed results`,async()=>{
    const h=pageHarness(file);await h.ctx.calculate3Omega();
    const data=h.data();assert.equal(data.length,401);
    const first=data[0],J=context.thermalIntegral(.02,5e-6,1e-6,.0005),prefactor=1/(10*Math.PI);
    assert(Math.abs(first.deltaT_in_phase-prefactor*J.re)<1e-12);
    assert(Math.abs(first.deltaT_out_of_phase-prefactor*J.im)<1e-12);
    assert(data.every(p=>p.integration_relative_error<=2e-5));
    assert.match(h.node('simulation-status').textContent,/numerical error/);
    assert.equal(h.node('save-data-button').disabled,false);
    h.ctx.exportDataFile();assert.equal(h.blobs.length,1);
    const exported=await h.blobs[0].text();
    assert.match(exported,/v1\.5\.1/);assert.match(exported,/Integration_rel_error_est/);
    const rows=exported.split('ln(2w)')[1].trim().split('\n').slice(1);
    assert.equal(rows.length,401);
    assert(rows.every(row=>Number(row.trim().split(/\s+/)[6])<=2e-5));
    h.ctx.thermalIntegral=()=>{throw Error('Integration did not converge.');};
    await h.ctx.calculate3Omega();assert.equal(h.data().length,0);
    assert.equal(h.node('save-data-button').disabled,true);assert.equal(h.node('save-pdf-button').disabled,true);
    assert.match(h.node('simulation-status').textContent,/did not converge/);
    h.ctx.exportDataFile();assert.equal(h.blobs.length,1);
  });
}
test('both pages and package identify release 1.5.1',()=>{
  assert.equal(require('../package.json').version,'1.5.1');
  for(const file of ['index.html','infinite_substrate.html']){
    const html=fs.readFileSync(path.join(root,file),'utf8');
    assert.match(html,/<title>.*v1\.5\.1/);
    for(const m of html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/gi))new vm.Script(m[1]);
  }
});
