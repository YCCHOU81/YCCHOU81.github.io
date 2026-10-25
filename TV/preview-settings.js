/* Overview sampling follows measured branch costs. Arithmetic precision and
   independent root checks remain the numerical workers' responsibility. */
(() => {
  "use strict";
  const profile=(id,label,quality,depth,budgetMs,reason)=>Object.freeze({id,label,quality,depth,density:1,budgetMs,meshBudgetMs:id==="tv-large"?45000:budgetMs,tailRefinement:0,reason});
  function recommend(input={}){
    const family=input.family??input.caseConfig?.seed??"lame";
    if(family==="lame"){
      const n=1+(input.counts??[0,0,0]).reduce((sum,value)=>sum+Number(value||0),0);
      if(n<=3)return profile("lame-small","Balanced Lamé overview",8,5,12000,"Small Lamé orders can show a deeper tail while keeping the interior responsive.");
      if(n<=11)return profile("lame-medium","Quick Lamé overview",8,3,12000,"Cusp continuation costs more than the compact interior at this order; start with a shorter checked tail.");
      return profile("lame-large","Large-order Lamé preview",4,3,10000,"Large Lamé orders need more work per root. Use a small interior mesh and extend the tails manually when useful.");
    }
    const c=input.caseConfig??{},weights=c.weights??[],total=weights.reduce((sum,value)=>sum+Number(value||0),0),largest=Math.max(0,...weights.map(Number));
    const operation=input.operation??c.operation??"repeated",repeated=input.isRepeated??operation==="repeated";
    if(family!=="bouquet"&&!repeated)return profile("tv-precomputed","Precomputed interior + cusp preview",20,operation==="seed"?5:3,12000,"The verified interior is already available. Local decorations use a shorter fresh cusp trace.");
    if(family==="bouquet"){
      if(total<=12&&largest<=8)return profile("bouquet-small","Balanced three-loop overview",8,5,15000,"Marked graph continuation and stable Fourier evaluation use a short first tail. Extend one edge in Limits & Zoom when you want to test the conjecture more closely.");
      if(total<=24&&largest<=15)return profile("bouquet-medium","Quick three-loop overview",4,3,25000,"Larger decorated cores use fewer interior samples and a short first tail. The extra time allowance lets graph initialization finish; checked points are retained.");
      return profile("bouquet-large","Large three-loop preview",4,3,35000,"High counts need more time to initialize the marked graph. Start with a compact overview, then extend a single edge; large resonance can make asymptotic convergence slow.");
    }
    if(total<=12&&largest<=8)return profile(c.connector?"tv-connector":"tv-mixed",c.connector?"Connector overview":"Mixed-decoration overview",8,3,12000,"The interior remains quick at these weights; fresh TV cusp transport is slower, so the first tail is shorter.");
    return profile("tv-large","Mixed / large TV preview",4,3,35000,"Use a small interior mesh and a short first tail. The allowance includes marked-graph initialization; difficult graphs may retain partial traces. Continue one edge in Limits & Zoom without discarding its checked points.");
  }
  window.AtlasPreview=Object.freeze({recommend});
})();
