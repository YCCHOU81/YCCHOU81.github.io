// Reuse completed, checked vertices of the same labelled graph and grid.
// This skips finished work only: fresh vertices still select their marked
// anchor and pass the ordinary continuation and endpoint checks.
const charts=['raw','zero','one','half'];
const pair=v=>Array.isArray(v)&&v.length===2&&v.every(x=>x!==null&&x!==''&&Number.isFinite(Number(x)));
function checked(branch,family,part){
 if(branch?.status!=='ok'||!pair(branch.tauExact)||Number(branch.tauExact[1])<=0)return false;
 const precision=Number(branch.precisionDigits),plot=family==='bouquet'&&part==='result'&&precision===15,tolerance=family==='bouquet'?(plot?1e-9:1e-12):1e-10;
 if(!Number.isFinite(precision)||precision<15)return false;
 if(plot){
  if(!(branch.plotChecked&&branch.twoFreshModels&&branch.pivotHealthy)||!Number.isFinite(Number(branch.plotErrorTolerance))||Number(branch.plotErrorTolerance)>tolerance)return false;
  if(!Number.isFinite(Number(branch.plotCorrectionTarget))||Number(branch.plotCorrectionTarget)>3e-12)return false;
 }else if(family==='bouquet'){
  if(precision<80||!branch.endpointChecked||!branch.pivotHealthy||Number(branch.lowRequestedDigits)<18||Number(branch.highRequestedDigits)<28||!Number.isFinite(Number(branch.lowRequestedDigits))||!Number.isFinite(Number(branch.highRequestedDigits)))return false;
 }else{
  const error=Number(branch.crossPrecisionChartError);
  if(precision<50||branch.crossPrecisionChartError==null||!Number.isFinite(error)||error<0||error>(family==='tetra'?1e-18:1e-10))return false;
 }
 if(!charts.every(name=>pair(branch.views?.[name])&&branch.viewErrors?.[name]!=null&&Number.isFinite(Number(branch.viewErrors[name]))&&Number(branch.viewErrors[name])>=0&&Number(branch.viewErrors[name])<=tolerance))return false;
 const [x,y]=branch.tauExact.map(Number),d=x*x+y*y,e=(1-x)**2+y*y,h=(1-2*x)**2+4*y*y;
 const expected={raw:[x,y],zero:[1-x/d,y/d],one:[(1-x)/e,y/e],half:[(x-2*d)/h,y/h]};
 return charts.every(name=>expected[name].every((value,i)=>Number.isFinite(value)&&Math.abs(value-Number(branch.views[name][i]))<=1e-11*(1+Math.abs(value))));
}
export function reusableMeshSamples(request,points,resolution,inset){
 const retained=new Map(),artifact=request.reuseMesh,result=artifact?.result,family=request.caseConfig?.seed;
 if(!['double','tetra','bouquet'].includes(family)||!artifact||artifact.caseId!==request.caseId||
    JSON.stringify(artifact.caseConfig)!==JSON.stringify(request.caseConfig)||
    Number(result?.resolution)!==resolution||Number(result?.inset)!==inset||
    !Array.isArray(result?.points)||!Array.isArray(result.samples)||result.points.length!==result.samples.length||result.points.length>points.length)return retained;
 const ids=new Map(points.map((p,i)=>[p.r+'|'+p.s,i]));
 for(let k=0;k<result.points.length;k++){
  const point=result.points[k],sample=result.samples[k],id=ids.get(point?.r+'|'+point?.s);
  if(id===undefined||retained.has(id)||sample?.source?.[0]!==point.r||sample?.source?.[1]!==point.s||
     !checked(sample.base,family,'base')||!checked(sample.result,family,'result'))return new Map();
  retained.set(id,sample);
 }
 return retained;
}
