"use strict";

/* Face labels are placed in the actual embedded ribbon faces, rather than at
   fixed locations.  This does not alter the numerical graph or its marking. */
globalThis.BouquetLayout=(()=>{
  const middle=[142,95],centers=[[90,33],[222,95],[90,157]];
  function direction(i){const p=centers[i-1],n=Math.hypot(p[0]-middle[0],p[1]-middle[1]);return[(p[0]-middle[0])/n,(p[1]-middle[1])/n];}
  function facePolygons(core,edges,positions,index){
    const curves=new Map();
    for(const e of edges){
      const a=positions[index.get(e.start)],b=positions[index.get(e.end)],loop=e.id.match(/^l([123])/);let points;
      if(loop){
        const u=direction(Number(loop[1])),v=[-u[1],u[0]],c1=[a[0]+78*u[0]-37*v[0],a[1]+78*u[1]-37*v[1]],c2=[b[0]+78*u[0]+37*v[0],b[1]+78*u[1]+37*v[1]];
        points=Array.from({length:49},(_,i)=>{const t=i/48,s=1-t;return[0,1].map(k=>s*s*s*a[k]+3*s*s*t*c1[k]+3*s*t*t*c2[k]+t*t*t*b[k]);});
      }else points=[a,b];
      curves.set(e.start,points);curves.set(e.end,points.slice().reverse());
    }
    return core.faces.map(f=>f.filter(d=>curves.has(d)).flatMap(d=>curves.get(d).slice(0,-1)));
  }
  function inside(p,poly){
    let yes=false;
    for(let i=0,j=poly.length-1;i<poly.length;j=i++){
      const a=poly[i],b=poly[j];
      if((a[1]>p[1])!==(b[1]>p[1])&&p[0]<(b[0]-a[0])*(p[1]-a[1])/(b[1]-a[1])+a[0])yes=!yes;
    }
    return yes;
  }
  function distanceSquared(p,poly){
    let out=Infinity;
    for(let i=0;i<poly.length;i++){
      const a=poly[i],b=poly[(i+1)%poly.length],dx=b[0]-a[0],dy=b[1]-a[1],norm=dx*dx+dy*dy;
      const t=norm?Math.max(0,Math.min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/norm)):0;
      out=Math.min(out,(p[0]-a[0]-t*dx)**2+(p[1]-a[1]-t*dy)**2);
    }
    return out;
  }
  function interiorLabel(poly){
    const xs=poly.map(p=>p[0]),ys=poly.map(p=>p[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
    let best=null,clearance=-1;
    const consider=p=>{if(inside(p,poly)){const d=distanceSquared(p,poly);if(d>clearance){best=p;clearance=d;}}};
    for(let i=0;i<poly.length;i++){
      const a=poly[i],b=poly[(i+1)%poly.length],dx=b[0]-a[0],dy=b[1]-a[1],n=Math.hypot(dx,dy);
      if(n)for(const offset of[-2,-.5,.5,2])consider([(a[0]+b[0])/2-offset*dy/n,(a[1]+b[1])/2+offset*dx/n]);
    }
    for(let i=0;i<=16;i++)for(let j=0;j<=16;j++)consider([xmin+(xmax-xmin)*(i+.5)/17,ymin+(ymax-ymin)*(j+.5)/17]);
    if(!best)throw new Error("The three-loop diagram has an unresolved face label");
    let radius=Math.max(xmax-xmin,ymax-ymin)/17;
    for(let round=0;round<5;round++){
      const center=best.slice();
      for(let i=-3;i<=3;i++)for(let j=-3;j<=3;j++)consider([center[0]+i*radius/3,center[1]+j*radius/3]);
      radius/=3;
    }
    return best;
  }
  function facePositions(core,edges,positions,index){
    const polygons=facePolygons(core,edges,positions,index),out={0:[260,20]};
    for(let i=1;i<4;i++)out[i]=interiorLabel(polygons[i]);
    return out;
  }
  return{facePositions,facePolygons,inside,interiorLabel,distanceSquared};
})();
