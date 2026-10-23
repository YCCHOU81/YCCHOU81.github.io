"use strict";

window.AtlasUI = (() => {
  const escape = value => String(value??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  function explainStop(result={}) {
    const raw=String(result.message??result.stopReason??""), stop=result.stop??{};
    const text=(String(stop.kind??"")+" "+String(result.stopReason??"")+" "+raw).toLowerCase();
    if(/exceeds.*(?:hemisphere|weight|local numerical budget)|unsupported|invalid.*(?:count|configuration)|outside.*(?:count|numerical budget)/.test(text))return {kind:"unsupported",title:"This configuration exceeds the numerical limits",explanation:"The count or face-weight limits were exceeded. This is a computational restriction, not a time limit or a mathematical obstruction.",advice:"Reduce the counts to the supported range shown beside the graph controls."};
    if(/time_budget|budget|deadline/.test(text)) return {kind:"time_budget",title:"Run time reached — samples kept",explanation:"This run reached its time allowance. The curve ends at a checked sample; it is not the mathematical endpoint.",advice:stop.resumeAvailable||result.resumeAvailable?"Continue the calculation to reuse its checkpoint. You can also choose fewer displayed samples for a new run.":"Try the shorter preview, fewer displayed samples, or a starting point closer to the target."};
    if(/cancelled|canceled|stopped by/.test(text)) return {kind:"cancelled",title:"Calculation cancelled — samples kept",explanation:"You stopped the calculation. The resolved points remain available for viewing and export.",advice:"Compute again to start a fresh run with the current settings."};
    if(/precision|consistency|relative.*check|branch.*reject|sector|independent.*agree|verification|overlap/.test(text)) return {kind:"verification",title:"A numerical check stopped the approach",explanation:"The next part did not pass a precision or marked-branch check. No unchecked continuation is drawn.",advice:"Try a shorter approach or move the starting point closer to the selected edge. Extra depth alone does not repair a failed check."};
    if(/unavailable|missing|unsupported|invalid|range|count/.test(text)) return {kind:"unsupported",title:"This calculation could not start",explanation:"The requested configuration or numerical data is not supported by this calculation.",advice:"Check the count limits and the detailed message below."};
    return {kind:"continuation",title:"The next sample could not be resolved",explanation:"The solver could not continue reliably. Earlier checked samples remain visible, with gaps where no point was resolved.",advice:"Try the shorter preview or move the source start closer to the target. A numerical stop does not imply that the surface fails to exist."};
  }
  function pageHref(page,search=location.search) {
    const params=new URLSearchParams(search);
    if(page==="explore.html" && params.has("fromR") && params.has("fromS")){
      params.set("r",params.get("fromR"));params.set("s",params.get("fromS"));
    }
    if(page==="limits.html" && params.has("r") && params.has("s") && !params.has("fromR")){
      params.set("fromR",params.get("r"));params.set("fromS",params.get("s"));
    }
    // Page navigation preserves the path without unexpectedly starting a solve.
    if(page==="limits.html")params.set("autocompute","0");
    return page+(params.size?"?"+params:"");
  }
  function refreshNavigation(search=location.search){
    document.querySelectorAll("a[data-page]").forEach(link=>{link.href=pageHref(link.dataset.page,search);});
  }
  function progressTable(container,traces,edgeLabels,makeLink){
    if(!container)return;
    const list=traces??[], partial=list.some(t=>t&&!["ok","computing"].includes(t.status));
    container.hidden=!list.some(Boolean);
    const details=container.querySelector("details"),wasOpen=details?.open??false;
    container.innerHTML=`<details${partial||wasOpen?" open":""}><summary>Cusp series progress · ${list.filter(t=>t?.status==="ok").length}/6 complete</summary><div class="table-scroll"><table class="series-progress-table"><thead><tr><th>Source edge</th><th>Series</th><th>Last resolved ε</th><th>Status</th><th>Inspect</th></tr></thead><tbody>${list.map((trace,i)=>{
      if(!trace)return "";
      const node=trace.nodes?.findLast(n=>n.status==="ok"&&n.result?.status==="ok");
      const state=trace.status==="ok"?"Complete":trace.status==="computing"?"Computing…":explainStop(trace).title;
      const title=escape(trace.message??trace.stop?.kind??state);
      return `<tr><td>${escape(edgeLabels[Math.floor(i/2)])}</td><td>${i%2+1}</td><td>${escape(node?.epsilon??"—")}</td><td title="${title}">${escape(state)}</td><td><a href="${escape(makeLink(trace,i))}">One-edge view →</a></td></tr>`;
    }).join("")}</tbody></table></div><p class="caption">Six series means two paths per cusp. A partial series stops at its last checked sample. Open a one-edge view for the reason and continuation options; Recompute runs the overview series again.</p></details>`;
  }
  function initWorkflowNavigation(){
    const nav=document.querySelector(".workflow-nav");
    if(!nav)return;
    const header=document.querySelector(".site-header"),main=nav.closest("main")||document.body;
    const steps=Array.from(nav.querySelectorAll('a[href^="#"]')).map(link=>{
      const heading=document.getElementById(link.getAttribute("href").slice(1));
      return heading?.classList.contains("workflow-step-heading")?{link,heading}:null;
    }).filter(Boolean);
    if(!steps.length)return;
    let frame=0,active=-1,headerHeight=-1,navHeight=-1;
    const schedule=()=>{if(!frame)frame=requestAnimationFrame(update);};
    function update(){
      frame=0;
      const nextHeader=Math.ceil(header?.getBoundingClientRect().height||0),nextNav=Math.ceil(nav.getBoundingClientRect().height);
      if(nextHeader!==headerHeight||nextNav!==navHeight){
        headerHeight=nextHeader;navHeight=nextNav;
        document.documentElement.style.setProperty("--workflow-header-height",headerHeight+"px");
        document.documentElement.style.setProperty("--workflow-nav-height",navHeight+"px");
        document.documentElement.style.setProperty("--workflow-scroll-offset",(headerHeight+navHeight+12)+"px");
      }
      // A step begins when its heading clears both sticky bars. Read only the
      // three or four cached headings; async graphs can move their positions.
      const threshold=headerHeight+navHeight+13;
      let next=0;
      for(let i=1;i<steps.length;i++)if(steps[i].heading.getBoundingClientRect().top<=threshold)next=i;
      const root=document.scrollingElement||document.documentElement,maxScroll=root.scrollHeight-window.innerHeight;
      if(maxScroll>1&&root.scrollTop>=maxScroll-2)next=steps.length-1;
      if(next===active)return;
      active=next;
      steps.forEach(({link},index)=>{if(index===active)link.setAttribute("aria-current","step");else link.removeAttribute("aria-current");});
      // Reveal the active tab on narrow screens without moving the page.
      const current=steps[active].link.getBoundingClientRect(),viewport=nav.getBoundingClientRect();
      if(current.left<viewport.left+8)nav.scrollLeft+=current.left-viewport.left-8;
      else if(current.right>viewport.right-8)nav.scrollLeft+=current.right-viewport.right+8;
    }
    window.addEventListener("scroll",schedule,{passive:true});
    window.addEventListener("resize",schedule,{passive:true});
    window.addEventListener("hashchange",schedule);
    window.addEventListener("pageshow",schedule);
    window.addEventListener("load",schedule,{once:true});
    window.addEventListener("atlas:cachedpage",schedule);
    if(typeof ResizeObserver!=="undefined"){
      const observer=new ResizeObserver(schedule);
      new Set([header,nav,main,document.body,...steps.map(step=>step.heading.parentElement)]).forEach(node=>{if(node)observer.observe(node);});
    }
    if(typeof MutationObserver!=="undefined"){
      // Decoration controls and numerical progress can change layout in place.
      // Ignore SVG geometry and our own aria-current/style updates.
      new MutationObserver(schedule).observe(main,{subtree:true,childList:true,characterData:true,attributes:true,attributeFilter:["hidden","open","class"]});
    }
    document.fonts?.ready.then(schedule);
    schedule();
  }
  document.addEventListener("DOMContentLoaded",()=>{
    refreshNavigation();
    initWorkflowNavigation();
    document.querySelectorAll("[data-copy-view]").forEach(button=>button.addEventListener("click",async()=>{
      const output=document.getElementById("view-action-status");
      try {await navigator.clipboard.writeText(location.href);if(output)output.textContent="Link copied with the current settings.";}
      catch {if(output)output.textContent="Copy the address from your browser to share these settings.";}
    }));
  });
  return {escape,explainStop,pageHref,refreshNavigation,progressTable};
})();
