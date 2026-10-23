// Research-only adaptive numerical precision check; not an interval certificate
// and not a substitute for continuation from a marked graph branch.
import{createEvaluator}from'./decimal-generic.mjs';
import{runCase}from'./decimal-case.mjs';

export function adaptiveConnector(spec,{horizontalDigits=12,verticalDigits=12,precisionIncrement=35,maxPrecision=900,maxAttempts=12}={}){
 const validation=createEvaluator(30);let epsilon;
 try{epsilon=new validation.D(spec.epsilon??'1e-12')}catch{return{success:false,precisionValidated:false,error:'epsilon must be a finite exact decimal string in (0,1]',stages:[],milliseconds:0}}
 if(!epsilon.isFinite()||epsilon.lte(0)||epsilon.gt(1))return{success:false,precisionValidated:false,error:'epsilon must be a finite exact decimal string in (0,1]',stages:[],milliseconds:0};
 const depth=epsilon.logarithm(10).neg().toNumber();
 maxPrecision=Math.min(970,Math.floor(maxPrecision));
 if(!Number.isInteger(precisionIncrement)||precisionIncrement<1||!Number.isInteger(maxAttempts)||maxAttempts<1||maxPrecision<25)return{success:false,precisionValidated:false,error:'Invalid precision budget',stages:[],milliseconds:0};
 let precision=spec.precision??Math.ceil(spec.edge===0?Math.max(75,3*depth+35):Math.max(80,depth+55));
 if(!Number.isInteger(precision)||precision<25||precision>maxPrecision)return{success:false,precisionValidated:false,error:'Required initial precision exceeds the configured budget',requiredPrecision:precision,maxPrecision,stages:[],milliseconds:0};
 let previous=null,guess=spec.guess,stages=[];
 for(let attempt=0;attempt<maxAttempts&&precision<=maxPrecision;attempt++,precision+=precisionIncrement){
  // Conservative attainable secant tolerance. A determinant with cancellation
  // may require more arithmetic digits; the next precision stage detects it.
  const correctionExponent=previous?new validation.D(previous.horizontalCorrection).e:0;
  let digits=previous?Math.max(30,horizontalDigits-correctionExponent+8):spec.edge===0?2*depth+horizontalDigits+6:35;
  if(!Number.isFinite(digits))digits=precision-depth-20;
  digits=Math.floor(Math.min(digits,precision-depth-22));
  const row=runCase({...spec,precision,digits,guess});stages.push(row);
  if(!row.success)continue;
  if(previous){
   const E=createEvaluator(Math.min(980,precision+20)),actual=E.c(row.w),old=E.c(previous.w),hc=new E.D(row.horizontalCorrection),ohc=new E.D(previous.horizontalCorrection);
   const relative=hc.minus(ohc).abs().div(hc.abs()),vertical=actual[1].minus(old[1]).abs();
   row.precisionComparison={horizontalRelativeError:relative.toString(),verticalAbsoluteError:vertical.toString()};
   if(relative.lt(new E.D(10).pow(-horizontalDigits))&&vertical.lt(new E.D(10).pow(-verticalDigits)))return{...row,precisionValidated:true,requestedHorizontalDigits:horizontalDigits,requestedVerticalDigits:verticalDigits,stages,milliseconds:stages.reduce((s,x)=>s+x.milliseconds,0)};
  }
  previous=row;guess=row.w;
 }
 return{...previous,success:false,precisionValidated:false,error:'Requested cross-precision accuracy was not obtained within the configured budget',maxPrecision,stages,milliseconds:stages.reduce((s,x)=>s+x.milliseconds,0)};
}
