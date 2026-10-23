// Checkpoints live in the numerical worker. Only a time-budget interruption is
// resumable: a rejected branch or an accuracy failure still needs inspection.
// The signature includes the labelled graph, exact source ray and display grid.
const checkpoints = new Map();
let serial = 0;
const canonical = value => Array.isArray(value) ? value.map(canonical) :
  value && typeof value === 'object' ? Object.fromEntries(Object.keys(value).sort().map(k => [k, canonical(value[k])])) : value;
// Compare exact decimal values without converting long source strings to a
// binary float. Numeric controls and their serialized decimal spellings agree.
const decimalKey = value => {
  const text=String(value),match=text.match(/^([+-]?)(\d*)(?:\.(\d*))?(?:e([+-]?\d+))?$/i);
  if(!match||!(match[2]||match[3]))return text;
  let digits=(match[2]+(match[3]??'')).replace(/^0+/,''),exponent=Number(match[4]??0)-(match[3]?.length??0);
  if(!digits)return '0';
  while(digits.endsWith('0')){digits=digits.slice(0,-1);exponent++;}
  return `${match[1]==='-'?'-':''}${digits}e${exponent}`;
};
export function traceSignature(request, family) {
  return JSON.stringify(canonical({family, caseConfig: request.caseConfig,
    path: request.path ?? {}, R: request.R ?? null, S: request.S ?? null, target: String(request.target ?? '1e-10'),
    startEpsilon: String(request.startEpsilon ?? 1),
    samplesPerDecade: Math.max(1, Math.min(5, Number(request.samplesPerDecade) || 2))}));
}
export function traceCheckpoint(request, family, create) {
  const signature = traceSignature(request, family);
  if (request.resumeToken || request.resume === true) {
    const saved = request.resumeToken ? checkpoints.get(String(request.resumeToken)) :
      [...checkpoints.values()].reverse().find(x => x.signature === signature && x.resumable);
    if (!saved) throw Error('This continuation checkpoint is no longer available. Compute the approach again.');
    if (saved.signature !== signature) throw Error('The case or path changed. Compute a new approach instead of resuming another branch.');
    if (!saved.resumable) throw Error('This trace stopped at a numerical verification failure and cannot be resumed automatically.');
    saved.attempts++;
    return {saved, resumed: true};
  }
  const token = `${family}-${Date.now().toString(36)}-${++serial}`;
  const saved = {token, signature, state: create(), resumable: false, attempts: 1};
  checkpoints.set(token, saved);
  while (checkpoints.size > 4) checkpoints.delete(checkpoints.keys().next().value);
  return {saved, resumed: false};
}
export function traceBudget(request) {
  return Math.max(250, Math.min(45000, Number(request.budgetMs) || 45000));
}
export function traceStop(saved, error, failedEpsilon, completed, total) {
  const message = String(error?.message ?? error), budget = /budget|deadline|time limit/i.test(message);
  saved.resumable = budget;
  return {kind: budget ? 'time_budget' : /precision|agree|rejected|match|sector|branch/i.test(message) ? 'verification' : 'numerical',
    message, failedEpsilon: failedEpsilon == null ? null : String(failedEpsilon),
    completed, total, retryable: budget, resumeAvailable: budget};
}
export function traceProgress(saved, resumed, stop = null) {
  return {stop, resumeToken: saved.resumable ? saved.token : null,
    resumed, attempts: saved.attempts, resumeAvailable: !!saved.resumable,
    ...(saved.state.extension ? {extended: true, reusedCount: saved.state.extension.reusedCount,
      newCount: Math.max(0, saved.state.nodes.filter(n => n.status === 'ok').length - saved.state.extension.reusedCount)} : {})};
}

// A serialized checked trace can seed a deeper request after page navigation or
// worker disposal. It is not a checkpoint: never import a live resume token.
// In particular, sparse local cusp roots do not establish marked continuation.
export function traceExtension(request, family) {
  if (!request.extendFrom || request.resumeToken || request.resume === true) return null;
  const previous = request.extendFrom.request, result = request.extendFrom.result;
  if (!previous || !result) throw Error('The saved trace has no complete request. Compute a new approach.');
  const identity = r => {
    if(r.seed&&r.caseConfig?.seed&&r.seed!==r.caseConfig.seed)throw Error('The saved seed label changed. Compute a new approach.');
    const seed=r.caseConfig?.seed??r.seed??family,path=r.path??{},start=path.start??[r.R??(seed==='bouquet'?'.29':'.58'),r.S??(seed==='double'?'.31':seed==='bouquet'?'.31':'.62')];
    return JSON.stringify(canonical({family,seed,
      caseId:r.caseId ?? r.caseConfig?.id ?? null,caseConfig:r.caseConfig,
      path:{kind:path.kind??'edge',index:Number(path.index??0),position:decimalKey(path.position??'.5'),start:start.map(decimalKey)},
      startEpsilon:decimalKey(r.startEpsilon??1),samplesPerDecade:Math.max(1,Math.min(5,Number(r.samplesPerDecade)||2))}));
  };
  if (identity(previous) !== identity(request)) throw Error('The graph, source path, starting epsilon or sampling density changed. Compute a new approach.');
  const oldTarget = Number(previous.target ?? '1e-10'), target = Number(request.target ?? '1e-10');
  if (!(oldTarget > target && target > 0)) throw Error('Trace extension requires a smaller final epsilon than the saved request.');
  if (['verification','numerical'].includes(result.stop?.kind)) throw Error('This trace stopped at a numerical or branch verification failure and cannot be extended automatically.');
  if (!['ok','partial','stopped'].includes(result.status)) throw Error('This trace has no safely checked prefix to extend.');
  const userStopped = result.status === 'stopped' && !result.stop && !/fail|reject|mismatch|verification|did not agree/i.test(result.message ?? '');
  if (result.status !== 'ok' && !['time_budget','cancelled','canceled','page_change','family_change'].includes(result.stop?.kind) && !result.cancelled && !userStopped)
    throw Error('This trace stopped at a numerical or branch verification failure and cannot be extended automatically.');
  if (family === 'bouquet' && result.tailBranchValidationStatus === 'pending')
    throw Error('The saved cusp tail has not matched the marked graph branch. Fill its gaps before increasing depth.');
  const nodes = JSON.parse(JSON.stringify(result.nodes ?? []));
  let reusedCount = 0, gap = false, lastEpsilon = Infinity;
  for (const node of nodes) {
    const epsilon = Number(node.epsilon);
    if (!(epsilon > 0 && epsilon < lastEpsilon)) throw Error('The saved sampling order is incomplete. Compute a new approach.');
    lastEpsilon = epsilon;
    if (node.status !== 'ok') { gap = true; continue; }
    if (gap) throw Error('The saved trace has a separated tail beyond an unresolved gap. Fill its gaps before increasing depth.');
    for (const part of ['base','result']) {
      const branch = node[part], tau = branch?.tauExact;
      if (branch?.status !== 'ok' || !Array.isArray(tau) || tau.length !== 2 || !tau.every(x => Number.isFinite(Number(x))) || !(Number(tau[1]) > 0))
        throw Error('The saved trace is missing exact checked period values. Compute a new approach.');
    }
    if (![node.rExact ?? node.r,node.sExact ?? node.s].every(x => x != null && Number.isFinite(Number(x))))
      throw Error('The saved trace is missing its exact source point. Compute a new approach.');
    reusedCount++;
  }
  if (!reusedCount) throw Error('The saved trace has no checked prefix to extend.');
  if (family === 'bouquet' && nodes.some(n => n.status === 'ok' && n.result?.engine?.includes('cusp')) && result.tailBranchValidationStatus !== 'matched')
    throw Error('The saved cusp tail has not matched the marked graph branch. Fill its gaps before increasing depth.');
  if (result.status === 'ok' && gap) throw Error('The saved complete trace contains an unresolved gap. Compute a new approach.');
  if (result.status === 'ok' && Number(nodes.at(-1).epsilon) !== oldTarget) throw Error('The saved complete trace does not reach its original final epsilon.');
  if (Number(nodes.at(-1).epsilon) < target || Number(nodes[0].epsilon) < oldTarget)
    throw Error('The saved trace does not belong to the requested epsilon interval.');
  return {nodes, reusedCount, last:nodes[reusedCount-1], result,
    oldTarget:String(previous.target ?? '1e-10')};
}

// Preserve every old display sample and its exact decimal spelling. New steps
// begin strictly after the old endpoint, rather than redistributing its grid.
export function extendTraceGrid(extension, target, density, power) {
  const wanted = extension.nodes.map(n => String(n.epsilon)), from = -Math.log10(Number(wanted.at(-1))), to = -Math.log10(Number(target));
  const count = Math.max(1, Math.ceil((to-from)*density));
  for (let k=1;k<=count;k++) wanted.push(k===count ? String(target) : power(-from-(to-from)*k/count));
  return wanted;
}

