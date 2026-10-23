/* Bounded, data-only browser cache. IndexedDB carries computed data between
   pages; a memory-only fallback keeps the site usable if storage is unavailable.
   The budget counts serialized UTF-16 strings, not total browser RAM. */
(() => {
  "use strict";
  const limitBytes = 16 * 1024 * 1024, maxEntries = 32;
  const entries = new Map();
  let bytes = 0, evictions = 0, db = null, storage = "memory", revision=0;
  let ready = null, pending = Promise.resolve();
  const id = (namespace, key) => JSON.stringify([String(namespace), String(key)]);
  const dirty = new Set();
  let hydrating = false;
  function remove(entryId) {
    const old = entries.get(entryId);
    if (!old) return false;
    bytes -= old.bytes; entries.delete(entryId); return true;
  }
  function failStorage() {
    storage = "memory";
    try { db?.close(); } catch (_) {}
    db = null;
  }
  function persist(changes = []) {
    if (!ready && !db) return;
    const order = [...entries.keys()];
    const localSizes=Object.fromEntries([...entries].map(([entryId,entry])=>[entryId,entry.bytes]));
    pending = pending.catch(() => {}).then(() => ready).then(() => new Promise(resolve => {
      if (!db) return resolve();
      try {
        const tx = db.transaction("entries", "readwrite"), store = tx.objectStore("entries");
        // Merge the small global index transactionally: open tabs share this
        // disk quota even when another tab has written since hydration.
        const index=store.get("@order");
        index.onsuccess=()=>{
          const previous=index.result??{},sizes={...(previous.sizes??{})};
          let ids=Array.isArray(previous.ids)?[...new Set(previous.ids)]:[];
          for(const entryId of ids)if(!Number.isFinite(sizes[entryId])&&Number.isFinite(localSizes[entryId]))sizes[entryId]=localSizes[entryId];
          const localOrder=order.filter(entryId=>ids.includes(entryId)),touched=new Set(localOrder);
          ids=ids.filter(entryId=>!touched.has(entryId)).concat(localOrder);
          for(const [entryId,entry]of changes){
            ids=ids.filter(value=>value!==entryId);
            if(entry){store.put({id:entryId,...entry});sizes[entryId]=entry.bytes;ids.push(entryId);}
            else{store.delete(entryId);delete sizes[entryId];}
          }
          ids=ids.filter(entryId=>{if(Number.isFinite(sizes[entryId])&&sizes[entryId]>=0)return true;store.delete(entryId);delete sizes[entryId];return false;});
          let total=ids.reduce((sum,entryId)=>sum+sizes[entryId],0);
          while(ids.length>maxEntries||total>limitBytes){const entryId=ids.shift();total-=sizes[entryId];delete sizes[entryId];store.delete(entryId);}
          store.put({id:"@order",ids,sizes});
        };
        tx.oncomplete = () => resolve();
        tx.onerror = tx.onabort = () => {failStorage(); resolve();};
      } catch (_) {failStorage(); resolve();}
    }));
  }
  function key(value) {
    const canonical = item => Array.isArray(item) ? item.map(canonical)
      : item && typeof item === "object"
        ? Object.fromEntries(Object.keys(item).sort().map(name => [name, canonical(item[name])])) : item;
    return JSON.stringify(canonical(value));
  }
  function forget(entryId) {
    revision++;
    if(hydrating)dirty.add(entryId); const removed = remove(entryId);
    persist([[entryId, null]]); return removed;
  }
  function loadRecords(rows, replace=false) {
    if(replace){entries.clear();bytes=0;}
    const savedOrder=rows.find(row=>row?.id==="@order")?.ids,order=Array.isArray(savedOrder)?savedOrder:[],ranks=new Map(order.map((entryId,index)=>[entryId,index]));
    const invalid=[];
    rows.filter(row=>row?.id!=="@order").sort((a,b)=>(ranks.get(a?.id)??-1)-(ranks.get(b?.id)??-1)).forEach(row=>{
      if(!row||typeof row.id!=="string"||typeof row.json!=="string"||typeof row.namespace!=="string"||typeof row.key!=="string"||row.id!==id(row.namespace,row.key)){if(row?.id)invalid.push([row.id,null]);return;}
      const size=2*(row.id.length+row.json.length);if(dirty.has(row.id))return;
      if(size>limitBytes){invalid.push([row.id,null]);return;}
      try{JSON.parse(row.json);}catch(_){invalid.push([row.id,null]);return;}
      remove(row.id);entries.set(row.id,{namespace:row.namespace,key:row.key,json:row.json,bytes:size});bytes+=size;
    });
    while(entries.size>maxEntries||bytes>limitBytes){const entryId=entries.keys().next().value;remove(entryId);invalid.push([entryId,null]);evictions++;}
    if(invalid.length)persist(invalid);
  }
  if (typeof indexedDB !== "undefined") {
    hydrating = true;
    ready = new Promise(resolve => {
      let settled = false;
      const done = () => {if (!settled) {settled = true; hydrating=false;dirty.clear();resolve();}};
      // A blocked database must not hold the controls hostage.
      const timeout = typeof setTimeout === "function" ? setTimeout(() => {failStorage(); done();}, 1500) : null;
      const finish = () => {if (timeout !== null) clearTimeout(timeout); done();};
      try {
        const open = indexedDB.open("lame-atlas-computed-v1", 1);
        open.onupgradeneeded = () => {
          if (!open.result.objectStoreNames.contains("entries")) open.result.createObjectStore("entries", {keyPath: "id"});
        };
        open.onerror = open.onblocked = () => {failStorage(); finish();};
        open.onsuccess = () => {
          if (settled) {open.result.close(); return;}
          db = open.result; storage = "indexeddb";
          db.onversionchange = () => failStorage();
          try {
            const tx = db.transaction("entries", "readonly"), request = tx.objectStore("entries").getAll();
            request.onsuccess = () => {
              try {
                loadRecords(request.result??[]);
              } catch (_) {failStorage();}
              finish();
            };
            request.onerror = tx.onerror = tx.onabort = () => {failStorage();finish();};
          } catch (_) {failStorage();finish();}
        };
      } catch (_) {failStorage();finish();}
    });
  }
  window.AtlasResults = Object.freeze({
    ready, key,
    set(namespace, cacheKey, value) {
      revision++;
      const entryId = id(namespace, cacheKey); if(hydrating)dirty.add(entryId);
      let json;
      try {json = JSON.stringify(value, (name,item) => name === "resumeToken" || name === "resumeAvailable" ? undefined : item);}
      catch (_) {forget(entryId);return false;}
      if(typeof json !== "string") {forget(entryId);return false;}
      const size=2*(entryId.length+json.length);remove(entryId);
      if(size > limitBytes) {persist([[entryId,null]]);return false;}
      const changes=[];
      while(entries.size >= maxEntries || bytes+size > limitBytes) {const entryId=entries.keys().next().value;remove(entryId);changes.push([entryId,null]);evictions++;}
      const entry={namespace:String(namespace),key:String(cacheKey),json,bytes:size};
      entries.set(entryId,entry);bytes+=size;changes.push([entryId,entry]);persist(changes);return true;
    },
    get(namespace,cacheKey) {
      const entryId=id(namespace,cacheKey),entry=entries.get(entryId);if(!entry)return null;
      entries.delete(entryId);entries.set(entryId,entry);
      try{return JSON.parse(entry.json);}catch(_){forget(entryId);return null;}
    },
    list(namespace) {
      const result=[];
      for(const entry of entries.values())if(entry.namespace===String(namespace))try{result.push({key:entry.key,value:JSON.parse(entry.json)});}catch(_){}
      return result;
    },
    drop(namespace,cacheKey){return forget(id(namespace,cacheKey));},
    clear(namespace) {
      revision++;
      const changes=[];
      for(const [entryId,entry] of entries)if(namespace===undefined||entry.namespace===String(namespace)){if(hydrating)dirty.add(entryId);remove(entryId);changes.push([entryId,null]);}
      persist(changes);
    },
    async flush(){await ready;persist();await pending;},
    async refresh(){
      await ready;await pending;if(!db)return;
      const before=revision;
      await new Promise(resolve=>{
        try{
          const tx=db.transaction("entries","readonly"),request=tx.objectStore("entries").getAll();
          request.onsuccess=()=>{try{if(revision===before)loadRecords(request.result??[],true);}catch(_){failStorage();}resolve();};
          request.onerror=tx.onerror=tx.onabort=()=>{failStorage();resolve();};
        }catch(_){failStorage();resolve();}
      });
    },
    stats(){return{entries:entries.size,bytes,limitBytes,maxEntries,evictions,storage};}
  });
  // Write the current view before ordinary atlas navigation, then hydrate the
  // destination cache before its numerical controller starts.
  if(typeof document !== "undefined" && document.addEventListener) {
    document.addEventListener("click",event=>{
      if(event.defaultPrevented||event.button>0||event.ctrlKey||event.metaKey||event.shiftKey||event.altKey)return;
      const link=event.target?.closest?.("a[href]");
      if(!link||link.target&&link.target!=="_self"||link.hasAttribute("download"))return;
      let next;try{next=new URL(link.href,location.href);}catch(_){return;}
      if(next.origin!==location.origin||next.pathname===location.pathname||!/(?:explore|limits|Lame|TV|mathematics)\.html$|\/overlap\/index\.html$/.test(next.pathname))return;
      event.preventDefault();
      window.dispatchEvent(new CustomEvent("atlas:beforepagechange"));
      window.AtlasResults.flush().finally(()=>{location.href=next.href;});
    });
    window.addEventListener?.("pagehide",()=>{window.dispatchEvent(new CustomEvent("atlas:beforepagechange"));window.AtlasResults.flush();});
    window.addEventListener?.("pageshow",event=>{if(event.persisted)window.AtlasResults.refresh().finally(()=>window.dispatchEvent(new CustomEvent("atlas:cachedpage")));});
  }
})();
