"use strict";
globalThis.GraphArcView = (() => {
  const NS="http://www.w3.org/2000/svg";
  const escape=x=>String(x).replaceAll('&','&amp;').replaceAll('<','&lt;').replaceAll('"','&quot;');
  const subs=n=>String(n).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[Number(d)]);
  const math=body=>`<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">${body}</math>`;
  const name=r=>r.kind==='axis'?`${r.reversed?'−':''}L${subs(r.i)}`:r.kind==='pair'?`${r.reversed?'−':''}2L${subs(r.i)}`:r.kind==='period'?`Period ${r.cycle==='b'?'τ':'1'} rotation`:`Axis h${subs(r.i)}`;
  function create({onClear=()=>{},onEdge=()=>{},onArc=()=>{}}={}){
    const get=id=>document.getElementById(id);
    get('arc-clear').addEventListener('click',onClear);
    let links=[];
    get('arc-edge-chips').addEventListener('click',event=>{
      const button=event.target.closest('[data-linked-arc]');
      const link=button&&links.find(item=>item.id===button.dataset.linkedArc);
      if(link)onArc(link.request);
    });
    for(const id of ['arc-primal','arc-dual']){
      const select=event=>{
        if(event.type==='keydown'&&!['Enter',' '].includes(event.key))return;
        const edge=event.target.closest('.graph-edge[data-edge]');
        if(!edge)return;
        event.preventDefault();onEdge(edge.dataset.edge);
      };
      get(id).addEventListener('click',select);get(id).addEventListener('keydown',select);
    }
    function update({result,scene,config,request,edgeLink=null}){
      document.querySelectorAll('.graph-view .is-arc-edge').forEach(e=>e.classList.remove('is-arc-edge','is-arc-negative'));
      document.querySelectorAll('.graph-view .is-reverse-edge').forEach(e=>e.classList.remove('is-reverse-edge'));
      document.querySelectorAll('.arc-edge-badge').forEach(e=>e.remove());
      const has=!!result||!!edgeLink;
      get('arc-linked-graphs').hidden=!has;
      get('arc-show-graphs').hidden=!has;
      get('arc-clear').hidden=!has;
      get('arc-selection-status').classList.toggle('has-selection',has);
      if(!has){links=[];get('arc-selection-title').textContent='Click a red arc or an edge to follow its length.';get('arc-selection-note').textContent='Arcs and concrete gluing edges are linked in both directions.';return;}
      let terms;
      if(edgeLink){
        const {edge,matches}=edgeLink;links=matches;
        const count=matches.length;
        get('arc-selection-title').textContent=`Selected edge ${edge.label} · ${count} matching arc sum${count===1?'':'s'}`;
        get('arc-selection-note').textContent=count?'Highlighted arcs use this concrete edge. Their badges show its signed coefficient.':`No displayed arc sum uses ${edge.label} in this chart continuation. Try another period or half-turn step.`;
        get('arc-clear').textContent='Clear edge';
        get('arc-linked-title').textContent=`${edge.label} → its arc sums`;
        get('arc-edge-equation').innerHTML=count?matches.map(link=>math(`${link.result.symbolMath}<mo>=</mo>${link.result.expressionMath}<mo>=</mo>${link.result.metricExpressionMath}`)).join(''):'<p class="arc-no-match">This edge has no contribution to the currently displayed arc sums.</p>';
        get('arc-edge-chips').innerHTML=matches.map(link=>`<button type="button" class="arc-edge-chip edge-arc-link${link.coefficient<0?' is-negative':''}" data-linked-arc="${escape(link.id)}" aria-label="Select ${escape(link.label)} and its full edge sum">${escape(link.label)} <span>${link.coefficient<0?'−':'+'}${Math.abs(link.coefficient)} × ${escape(edge.label)}</span></button>`).join('');
        get('arc-edge-note').textContent=`${edge.label} is emphasized in both drawings. ${matches.some(link=>link.request.kind==='axis-location')?'A single half-turn link locates its axis; its jump also depends on the incoming point.':'Select a listed arc to see its full edge sum. Matching uses concrete edges, even when lengths coincide.'}`;
        terms=[{edgeId:edge.id,label:edge.label,coefficient:1}];
      }else{
        links=[];
        get('arc-clear').textContent='Clear arc';
        get('arc-selection-title').textContent=`Selected: ${name(request)} · ${result.terms.length} contributing edge${result.terms.length===1?'':'s'}`;
        get('arc-selection-note').textContent=request.kind==='axis-location'?(result.terms.length?'The highlight locates the half-turn axis; the jump also depends on the incoming point.':'This axis is the chosen angular origin; its jump still depends on the incoming point.'):'Solid edges add; dashed edges subtract. A coefficient counts repeated use.';
        get('arc-linked-title').textContent=`Edges for ${name(request)}`;
        get('arc-edge-equation').innerHTML=math(`${result.symbolMath}<mo>=</mo>${result.expressionMath}<mo>=</mo>${result.metricExpressionMath}`);
        get('arc-edge-chips').innerHTML=result.terms.length?result.terms.map(t=>`<span class="arc-edge-chip${t.coefficient<0?' is-negative':''}">${t.coefficient<0?'−':'+'}${Math.abs(t.coefficient)===1?'':`${Math.abs(t.coefficient)} × `}${escape(t.label)}</span>`).join(''):'<span>The signed edge contributions cancel.</span>';
        get('arc-edge-note').textContent=request.kind==='axis-location'?result.note:`ℓⱼ is the length of Eⱼ. Solid edges add; dashed edges subtract. ${['period','pair'].includes(request.kind)?'The angle is unreduced; full turns do not change its endpoint.':''} Click an edge in either drawing to find its arc sums.`;
        terms=result.terms;
      }
      const key=config.id;
      if(get('arc-primal').dataset.graphKey!==key){
        get('arc-primal').innerHTML=scene.primalSvg;get('arc-dual').innerHTML=scene.dualSvg;
        for(const id of ['arc-primal','arc-dual']){
          const view=get(id);view.dataset.graphKey=key;
          view.querySelector('svg').setAttribute('role','group');
          view.querySelectorAll('[data-edge],[data-hemisphere],[data-pole]').forEach(e=>{
            const enabled=e.matches('.graph-edge[data-edge]');e.dataset.interactive=String(enabled);
            e.setAttribute('tabindex',enabled?'0':'-1');
            if(enabled){const item=scene.metadata.edges.find(item=>item.id===e.dataset.edge);e.setAttribute('role','button');e.setAttribute('aria-label',`Select ${item.label} and highlight its arc sums`);}
            else{e.removeAttribute('role');e.removeAttribute('aria-pressed');}
          });
        }
      }
      for(const t of terms){
        document.querySelectorAll(`.graph-view .graph-edge[data-edge="${t.edgeId}"]`).forEach(group=>{
          group.classList.add('is-arc-edge');group.classList.toggle('is-arc-negative',t.coefficient<0);
          group.classList.toggle('is-reverse-edge',!!edgeLink);
          const path=group.querySelector('.edge-line');
          if(!path)return;
          const p=path.getPointAtLength(path.getTotalLength()/2),label=document.createElementNS(NS,'text');
          label.setAttribute('x',p.x);label.setAttribute('y',p.y-8);label.setAttribute('class','arc-edge-badge');label.setAttribute('text-anchor','middle');
          label.textContent=edgeLink?t.label:`${t.coefficient<0?'−':'+'}${Math.abs(t.coefficient)===1?'':`${Math.abs(t.coefficient)}×`}${t.label}`;
          group.append(label);
        });
      }
    }
    return {update};
  }
  return {create};
})();
