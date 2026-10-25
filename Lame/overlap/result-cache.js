/* Bounded Overlap Atlas results. Keep database/store names and record keys stable:
 * existing checked interiors and tails must survive controller refactoring.
 * The main atlas has a separate cache and memory budget.
 */
(function (root) {
  "use strict";

  function create(options = {}) {
    const databaseName = options.databaseName ?? "lame-atlas-overlap-v1";
    const maxBytes = options.maxBytes ?? 8 * 1024 * 1024;
    const maxItems = options.maxItems ?? 600;
    const memory = new Map();
    let database = null;
    let bytes = 0;

    function remember(id, record, when = Date.now()) {
      const size = 2 * (id.length + JSON.stringify(record).length);
      const previous = memory.get(id);
      if (previous) bytes -= previous.size;
      memory.delete(id);
      memory.set(id, {record, size, when});
      bytes += size;

      const removed = [];
      while (bytes > maxBytes || memory.size > maxItems) {
        const first = memory.keys().next().value;
        bytes -= memory.get(first).size;
        memory.delete(first);
        removed.push(first);
      }
      return removed;
    }

    const ready = new Promise(resolve => {
      if (typeof root.indexedDB === "undefined") {
        resolve();
        return;
      }

      let settled = false;
      const finish = () => {
        if (!settled) {
          settled = true;
          resolve();
        }
      };
      const timer = root.setTimeout(finish, 1500);
      const finishLoading = () => {
        root.clearTimeout(timer);
        finish();
      };

      try {
        const request = root.indexedDB.open(databaseName, 1);
        request.onupgradeneeded = () => request.result.createObjectStore("meshes", {keyPath: "key"});
        request.onerror = request.onblocked = finishLoading;
        request.onsuccess = () => {
          database = request.result;
          const get = database.transaction("meshes", "readonly").objectStore("meshes").getAll();
          get.onsuccess = () => {
            for (const row of get.result.sort((a, b) => a.when - b.when)) {
              if (row?.record?.counts && Array.isArray(row.record.tau) && !memory.has(row.key)) {
                remember(row.key, row.record, row.when);
              }
            }

            const stale = get.result.filter(row => !memory.has(row.key));
            if (stale.length) {
              try {
                const store = database.transaction("meshes", "readwrite").objectStore("meshes");
                for (const row of stale) store.delete(row.key);
              } catch (_) {
                // Storage can be unavailable; checked in-memory results remain usable.
              }
            }
            finishLoading();
          };
          get.onerror = finishLoading;
        };
      } catch (_) {
        finishLoading();
      }
    });

    return {
      ready,
      get(id) {
        const row = memory.get(id);
        if (!row) return null;
        memory.delete(id);
        memory.set(id, row);
        return row.record;
      },
      put(id, record) {
        const removed = remember(id, record);
        if (!database) return;
        try {
          const store = database.transaction("meshes", "readwrite").objectStore("meshes");
          store.put({key: id, record, when: Date.now()});
          for (const old of removed) store.delete(old);
        } catch (_) {
          // A write failure must not interrupt drawing or numerical work.
        }
      },
      get persistent() {
        return !!database;
      }
    };
  }

  root.OverlapResultCache = Object.freeze({create});
})(typeof globalThis !== "undefined" ? globalThis : this);
