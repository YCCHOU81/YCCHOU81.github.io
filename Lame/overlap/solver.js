/* Lamé classic-worker entry. Shared implementation lives in ../workers/.
 * Preserve the standalone solve/compute protocol as well as Overlap imports. */
"use strict";
importScripts("../workers/decimal.js?v=20261025-maintenance1",
  "../workers/atlas-math.js?v=20261025-maintenance1",
  "../workers/precision-engine.js?v=20261025-maintenance1",
  "../workers/attachment-engine.js?v=20261025-maintenance1",
  "../workers/attachment-anchors.js?v=20261025-maintenance1");
const classicalPrecision=createClassicalPrecision({reuseAcceptedEvaluation:true});
const atlasAttachment=createAttachmentEngine(classicalPrecision);
atlasAttachment.setAnchorData(attachmentAnchors);
importScripts("../workers/interior-messages.js?v=20261025-maintenance1");
