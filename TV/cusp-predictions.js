"use strict";

/* Simulation annotations only. None of these functions changes a numerical
   root, inserts an endpoint, or feeds an asymptotic value into the solver. */
window.TVCuspPredictions = (() => {
  const gcd = (a,b) => { a=Math.abs(a);b=Math.abs(b);while(b){const t=a%b;a=b;b=t;}return a||1; };
  function rational(numerator,denominator) {
    if(!Number.isSafeInteger(numerator)||!Number.isSafeInteger(denominator)||denominator<=0)return null;
    const d=gcd(numerator,denominator);numerator/=d;denominator/=d;
    return {numerator,denominator,value:numerator/denominator,label:denominator===1?String(numerator):`${numerator}/${denominator}`};
  }
  function chartValues(node,branch,chart) {
    const part=node?.[branch];
    if(part?.status!=="ok"||(branch==="result"&&node.status!=="ok"))return null;
    const exact=part.views?.[chart];
    if(exact&&exact.length===2&&Number.isFinite(Number(exact[0]))&&Number(exact[1])>0)return exact;
    // The worker normally supplies the opened chart directly. Do not invert
    // tiny, ordinary-precision raw coordinates to manufacture a deep cusp.
    if(chart!=="raw")return null;
    const raw=part.tauExact??(part.tau?[part.tau.re,part.tau.im]:null);
    return raw&&Number.isFinite(Number(raw[0]))&&Number(raw[1])>0?raw:null;
  }
  function numericalCandidate(points,{weights=[]}={}) {
    // Use a consecutive resolved suffix, never jump across an unresolved gap.
    const tail=[];
    for(let i=points.length-1;i>=0&&tail.length<6;i--){
      const p=points[i];if(!p||!Number.isFinite(Number(p[0]))||!(Number(p[1])>0))break;
      tail.unshift([Number(p[0]),Number(p[1])]);
    }
    if(tail.length<4)return null;
    const recent=tail.slice(-4),xs=recent.map(p=>p[0]),ys=recent.map(p=>p[1]);
    if(ys.at(-1)-ys[0]<.15||ys.some((y,i)=>i&&y<ys[i-1]-.05))return null;
    const count=weights.reduce((sum,n)=>sum+Math.max(0,Number(n)||0),0);
    const maxDenominator=Math.min(4096,Math.max(64,16*(count+1)));
    const fits=[{xs,samples:recent.length,extrapolated:false}];
    // A slowly narrowing cusp can have a visible horizontal correction long
    // after its geometric decay rate has settled. Three independent Aitken
    // estimates may support a candidate sooner. They annotate the plot only;
    // no extrapolated point is ever added to the trace.
    if(tail.length>=6){
      const estimates=[],ratios=[];
      for(let i=2;i<tail.length;i++){
        const d1=tail[i-1][0]-tail[i-2][0],d2=tail[i][0]-tail[i-1][0],q=d2/d1;
        if(!Number.isFinite(q)||q<=0||q>=.96||Math.abs(d1-d2)<1e-16){estimates.length=0;break;}
        estimates.push(tail[i][0]+d2*q/(1-q));ratios.push(q);
      }
      if(estimates.length>=3&&Math.max(...ratios)-Math.min(...ratios)<.04&&Math.abs(estimates.at(-1)-xs.at(-1))<.01)
        fits.push({xs:estimates.slice(-3),samples:tail.length,extrapolated:true});
    }
    for(const fit of fits){
      const spread=Math.max(...fit.xs)-Math.min(...fit.xs);
      if(spread>(fit.extrapolated?1e-7:1e-6))continue;
      const tolerance=Math.max(2e-10,Math.min(1e-7,2*spread)),x=fit.xs.at(-1);
      // Choose the simplest stable rational within the data's resolution,
      // rather than always returning a high-denominator approximation.
      for(let d=1;d<=maxDenominator;d++){
        const n=Math.round(x*d),r=rational(n,d);
        if(!r||Math.abs(x-r.value)>tolerance||fit.xs.some(v=>Math.abs(v-r.value)>4*tolerance))continue;
        return {...r,status:"numerical-candidate",kind:"Numerical rational candidate",samples:fit.samples,maxDenominator,extrapolated:fit.extrapolated,
          uncertainty:Math.max(spread,Math.abs(x-r.value)),detail:`Inferred from ${fit.samples} consecutive checked samples${fit.extrapolated?", with three consistent geometric-tail extrapolations":""}; denominators up to ${maxDenominator} were tested. The reference is a numerical conjecture, not an exact-limit proof${fit.extrapolated?"; the drawn trace contains only computed roots":""}.`};
      }
    }
    return null;
  }
  function reference(info,points=[],options={}) {
    if(!info||options.chart&&info.chart!==options.chart)return null;
    if(info.known!==false&&Number.isFinite(info.x0)){
      let r=rational(info.numerator??info.xNumerator,info.denominator??info.xDenominator);
      if(!r&&typeof info.xLabel==="string"){
        const pair=info.xLabel.split("/").map(Number);r=rational(pair[0],pair[1]??1);
      }
      if(r)return {...r,status:"proposed",kind:info.countFormula?"Count-formula conjecture":"Proposed graph limit",info,
        detail:info.countFormula?`${info.countFormula.name}. ${info.countFormula.evidence}`:"Graph formula for this labelled decoration. Compare it with the computed trace."};
    }
    const inferred=numericalCandidate(points,options);
    return inferred?{...inferred,info}:null;
  }
  function fractionHTML(ref) {
    if(!ref)return "<span class=\"tv-rational-unresolved\">Not resolved yet</span>";
    return `<math aria-label="${ref.label}">${ref.denominator===1?`<mn>${ref.numerator}</mn>`:`<mfrac><mn>${ref.numerator}</mn><mn>${ref.denominator}</mn></mfrac>`}</math>`;
  }
  return {rational,chartValues,numericalCandidate,reference,fractionHTML};
})();
