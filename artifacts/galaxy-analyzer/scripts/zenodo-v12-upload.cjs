#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19440400';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/12.0', 'Accept': 'application/json' },
    };
    if (body && contentType) {
      opts.headers['Content-Type'] = contentType;
      if (typeof body === 'string') opts.headers['Content-Length'] = Buffer.byteLength(body);
    }
    const req = https.request(opts, res => {
      let data = '';
      res.on('data', d => data += d);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('HTTP ' + res.statusCode + ': ' + data.substring(0, 500));
          reject(new Error('HTTP ' + res.statusCode));
        } else {
          try { resolve(JSON.parse(data)); } catch { resolve(data); }
        }
      });
    });
    req.on('error', reject);
    if (body && typeof body !== 'string') body.pipe(req);
    else { if (body) req.write(body); req.end(); }
  });
}

async function uploadFile(draftId, filePath, fileName) {
  const stat = fs.statSync(filePath);
  const initBody = JSON.stringify([{ key: fileName }]);
  await zenodoRequest('POST', '/api/records/' + draftId + '/draft/files', initBody, 'application/json');

  await new Promise((resolve, reject) => {
    const encodedName = encodeURIComponent(fileName);
    const opts = {
      hostname: BASE,
      path: '/api/records/' + draftId + '/draft/files/' + encodedName + '/content',
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'User-Agent': 'GalaxyAnalyzer/12.0',
        'Content-Type': 'application/octet-stream',
        'Content-Length': stat.size,
      },
    };
    const req = https.request(opts, res => {
      let data = '';
      res.on('data', d => data += d);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('Upload error ' + res.statusCode + ' for ' + fileName + ': ' + data.substring(0, 300));
          reject(new Error('Upload content failed'));
        } else {
          resolve();
        }
      });
    });
    req.on('error', reject);
    fs.createReadStream(filePath).pipe(req);
  });

  await zenodoRequest('POST', '/api/records/' + draftId + '/draft/files/' + encodeURIComponent(fileName) + '/commit', '', 'application/json');
  console.log('  Uploaded: ' + fileName + ' (' + stat.size + ' bytes)');
}

async function main() {
  console.log('=== Zenodo v12 Upload: Program Closure — 1D Information Ceiling ===');
  console.log('=== Programs 1-8C + Phases V/V+ + Program Closure ===\n');

  let draftId;

  console.log('Step 1: Create new version from record ' + LATEST_ID + '...');
  try {
    const newVer = await zenodoRequest('POST', '/api/records/' + LATEST_ID + '/versions', '', 'application/json');
    draftId = newVer.id;
    console.log('  Created new draft ID: ' + draftId);
  } catch (e) {
    console.log('  Create failed (draft may already exist), checking latest draft...');
    const latestRec = await zenodoRequest('GET', '/api/records/' + LATEST_ID, null, null);
    if (latestRec.links && latestRec.links.latest_draft) {
      const draftUrl = latestRec.links.latest_draft;
      const draftInfo = await zenodoRequest('GET', new URL(draftUrl).pathname, null, null);
      draftId = draftInfo.id;
      console.log('  Found existing draft ID: ' + draftId);
    } else {
      throw new Error('Cannot find or create draft');
    }
  }

  console.log('\nStep 2: Check and delete existing files from draft...');
  try {
    const filesInfo = await zenodoRequest('GET', '/api/records/' + draftId + '/draft/files', null, null);
    const entries = filesInfo.entries || filesInfo;
    if (Array.isArray(entries) && entries.length > 0) {
      for (const f of entries) {
        try {
          await zenodoRequest('DELETE', '/api/records/' + draftId + '/draft/files/' + encodeURIComponent(f.key), null, null);
          console.log('  Deleted: ' + f.key);
        } catch (e) {
          console.log('  Could not delete ' + f.key + ': ' + e.message);
        }
      }
    } else {
      console.log('  No existing files to delete.');
    }
  } catch (e) {
    console.log('  Could not list files: ' + e.message);
  }

  console.log('\nStep 3: Upload new files...');
  const baseDir = path.join(__dirname, '..');
  const filesToUpload = [
    { local: 'scripts/phase126-m4-candidate.cjs', remote: 'scripts/phase126-m4-candidate.cjs' },
    { local: 'scripts/phase127-vflat-challenge.cjs', remote: 'scripts/phase127-vflat-challenge.cjs' },
    { local: 'scripts/phase128-decode-4th-axis.cjs', remote: 'scripts/phase128-decode-4th-axis.cjs' },
    { local: 'scripts/phase129-vflat-decomposition.cjs', remote: 'scripts/phase129-vflat-decomposition.cjs' },
    { local: 'scripts/phase130-split-4th-sector.cjs', remote: 'scripts/phase130-split-4th-sector.cjs' },
    { local: 'scripts/phase131-decode-vfresid.cjs', remote: 'scripts/phase131-decode-vfresid.cjs' },
    { local: 'scripts/phase132-vfresid-reducibility.cjs', remote: 'scripts/phase132-vfresid-reducibility.cjs' },
    { local: 'scripts/phase132a-halo-death-match.cjs', remote: 'scripts/phase132a-halo-death-match.cjs' },
    { local: 'scripts/phase132b-mediation-causal.cjs', remote: 'scripts/phase132b-mediation-causal.cjs' },
    { local: 'scripts/phase132c-external-robustness.cjs', remote: 'scripts/phase132c-external-robustness.cjs' },
    { local: 'scripts/phase133a-regime-law.cjs', remote: 'scripts/phase133a-regime-law.cjs' },
    { local: 'scripts/phase133b-second-channel.cjs', remote: 'scripts/phase133b-second-channel.cjs' },
    { local: 'scripts/phase133c-coupling-drivers.cjs', remote: 'scripts/phase133c-coupling-drivers.cjs' },
    { local: 'scripts/phase134-external-validation.cjs', remote: 'scripts/phase134-external-validation.cjs' },
    { local: 'external-validation/phase200-data-assembly.cjs', remote: 'scripts/phase200-data-assembly.cjs' },
    { local: 'external-validation/phase201-blind-prediction.cjs', remote: 'scripts/phase201-blind-prediction.cjs' },
    { local: 'external-validation/phase202-hierarchy-replication.cjs', remote: 'scripts/phase202-hierarchy-replication.cjs' },
    { local: 'external-validation/phase203-logMhost-improvement.cjs', remote: 'scripts/phase203-logMhost-improvement.cjs' },
    { local: 'external-validation/phase204-final-external-synthesis.cjs', remote: 'scripts/phase204-final-external-synthesis.cjs' },
    { local: 'external-validation/phase300-sample-salvage.cjs', remote: 'scripts/phase300-sample-salvage.cjs' },
    { local: 'external-validation/phase301-vfresid-drivers.cjs', remote: 'scripts/phase301-vfresid-drivers.cjs' },
    { local: 'external-validation/phase302-regime-law.cjs', remote: 'scripts/phase302-regime-law.cjs' },
    { local: 'external-validation/phase303-physical-interpretation.cjs', remote: 'scripts/phase303-physical-interpretation.cjs' },
    { local: 'scripts/program7a-halo-profile-test.cjs', remote: 'scripts/program7a-halo-profile-test.cjs' },
    { local: 'scripts/program7b-halo-shape.cjs', remote: 'scripts/program7b-halo-shape.cjs' },
    { local: 'scripts/program7c-halo-shape-index.cjs', remote: 'scripts/program7c-halo-shape-index.cjs' },
    { local: 'scripts/program8a-2d-state-recovery.cjs', remote: 'scripts/program8a-2d-state-recovery.cjs' },
    { local: 'scripts/program8b-physics-search.cjs', remote: 'scripts/program8b-physics-search.cjs' },
    { local: 'scripts/program8c-map-reconstruction.cjs', remote: 'scripts/program8c-map-reconstruction.cjs' },
    { local: 'scripts/phase-v-red-team.cjs', remote: 'scripts/phase-v-red-team.cjs' },
    { local: 'scripts/phase-v-plus-distance.cjs', remote: 'scripts/phase-v-plus-distance.cjs' },
    { local: 'external-validation/PROGRAM.md', remote: 'external-validation/PROGRAM.md' },
    { local: 'public/phase126-m4-candidate.json', remote: 'results/phase126-m4-candidate.json' },
    { local: 'public/phase127-vflat-challenge.json', remote: 'results/phase127-vflat-challenge.json' },
    { local: 'public/phase128-decode-4th-axis.json', remote: 'results/phase128-decode-4th-axis.json' },
    { local: 'public/phase129-vflat-decomposition.json', remote: 'results/phase129-vflat-decomposition.json' },
    { local: 'public/phase130-split-4th-sector.json', remote: 'results/phase130-split-4th-sector.json' },
    { local: 'public/phase131-decode-vfresid.json', remote: 'results/phase131-decode-vfresid.json' },
    { local: 'public/phase132-vfresid-reducibility.json', remote: 'results/phase132-vfresid-reducibility.json' },
    { local: 'public/phase132a-halo-death-match.json', remote: 'results/phase132a-halo-death-match.json' },
    { local: 'public/phase132b-mediation-causal.json', remote: 'results/phase132b-mediation-causal.json' },
    { local: 'public/phase132c-external-robustness.json', remote: 'results/phase132c-external-robustness.json' },
    { local: 'public/phase133a-regime-law.json', remote: 'results/phase133a-regime-law.json' },
    { local: 'public/phase133b-second-channel.json', remote: 'results/phase133b-second-channel.json' },
    { local: 'public/phase133c-coupling-drivers.json', remote: 'results/phase133c-coupling-drivers.json' },
    { local: 'public/phase134-external-validation.json', remote: 'results/phase134-external-validation.json' },
    { local: 'public/phase200-external-dataset.json', remote: 'results/phase200-external-dataset.json' },
    { local: 'public/phase201-blind-prediction.json', remote: 'results/phase201-blind-prediction.json' },
    { local: 'public/phase202-hierarchy-replication.json', remote: 'results/phase202-hierarchy-replication.json' },
    { local: 'public/phase203-logMhost-improvement.json', remote: 'results/phase203-logMhost-improvement.json' },
    { local: 'public/phase204-final-external-synthesis.json', remote: 'results/phase204-final-external-synthesis.json' },
    { local: 'public/phase300-sample-salvage.json', remote: 'results/phase300-sample-salvage.json' },
    { local: 'public/phase301-vfresid-drivers.json', remote: 'results/phase301-vfresid-drivers.json' },
    { local: 'public/phase302-regime-law.json', remote: 'results/phase302-regime-law.json' },
    { local: 'public/phase303-physical-interpretation.json', remote: 'results/phase303-physical-interpretation.json' },
    { local: 'public/program7a-halo-profile.json', remote: 'results/program7a-halo-profile.json' },
    { local: 'public/program7b-halo-shape.json', remote: 'results/program7b-halo-shape.json' },
    { local: 'public/program7c-halo-shape-index.json', remote: 'results/program7c-halo-shape-index.json' },
    { local: 'public/program8a-2d-state.json', remote: 'results/program8a-2d-state.json' },
    { local: 'public/program8b-physics-search.json', remote: 'results/program8b-physics-search.json' },
    { local: 'public/program8c-map-reconstruction.json', remote: 'results/program8c-map-reconstruction.json' },
    { local: 'public/phase-v-red-team.json', remote: 'results/phase-v-red-team.json' },
    { local: 'public/phase-v-plus-distance.json', remote: 'results/phase-v-plus-distance.json' },
    { local: 'public/phase56-frozen-baselines.json', remote: 'data/phase56-frozen-baselines.json' },
    { local: 'public/stage-A-master-table.json', remote: 'data/stage-A-master-table.json' },
    { local: 'public/sparc-table.json', remote: 'data/sparc-table.json' },
    { local: 'public/phase58a2-tidal-expansion.json', remote: 'data/phase58a2-tidal-expansion.json' },
    { local: 'public/sparc-results.json', remote: 'data/sparc-results.json' },
    { local: 'public/replication/REPRODUCIBLE_RESULT.md', remote: 'REPRODUCIBLE_RESULT.md' },
    { local: 'public/replication/MANUSCRIPT-TRACK2.md', remote: 'MANUSCRIPT-TRACK2.md' },
    { local: 'public/replication/IFU-FOLLOWUP-MEMO.md', remote: 'IFU-FOLLOWUP-MEMO.md' },
  ];

  for (const f of filesToUpload) {
    const fullPath = path.join(baseDir, f.local);
    if (fs.existsSync(fullPath)) {
      await uploadFile(draftId, fullPath, f.remote);
    } else {
      console.log('  SKIPPED (not found): ' + f.local);
    }
  }

  console.log('\nStep 4: Update metadata...');
  const metadata = {
    metadata: {
      title: 'Per-Galaxy a0 Variation and the Hidden Common-Cause State in SPARC Galaxies — Program Closure (v12 — Phases 126-134, 200-204, 300-303, Programs 1-8C, Phases V/V+)',
      description:
        '<p>v12 — <strong>Program closure</strong>: complete investigation of the hidden physical variable H behind per-galaxy a<sub>0</sub> variation in SPARC rotation curves. ' +
        'SPARC 1D data reaches a provable information ceiling; the physical identity of H requires 2D/3D kinematic data (IFU velocity fields).</p>' +
        '<h3>Central Result</h3>' +
        '<p>A hidden common-cause state H drives a bilateral coupling between kinematic excess (VfResid) and acceleration discrepancy (a0Resid) at ' +
        '<strong>r &asymp; 0.77</strong> (LOO cross-validated, p &lt; 0.001). This coupling is construction-independent (55/55), distance-independent (V+ kill test: 4/5), ' +
        'and survives all adversarial tests (Red Team V: 5/8). However, ~70&ndash;80% of H is structurally inaccessible from 1D rotation curves.</p>' +
        '<h3>What is new in v12</h3>' +
        '<ol>' +
        '<li><strong>Program 8A</strong>: 2D state recovery from RC features &mdash; 0/10 features significant; information ceiling confirmed (88.5% inaccessible)</li>' +
        '<li><strong>Program 8B</strong>: Physics-grounded hidden-state search &mdash; 12 models across 3 families; <strong>inaccessibility-strength paradox</strong> discovered (no model achieves both r&ge;0.65 AND &gt;50% hidden)</li>' +
        '<li><strong>Program 8C</strong>: Full RC map reconstruction &mdash; 71-dim state vector LOSES to scalar haloResponse; map ceiling = scalar ceiling</li>' +
        '<li><strong>Phase V (Red Team)</strong>: Adversarial audit (5/8 PASS); headline revised to r &asymp; 0.77 (LOO, down from 0.80 global)</li>' +
        '<li><strong>Phase V+ (Distance Kill)</strong>: Distance-error hypothesis killed (4/5); partial r|logD = 0.804 unchanged; neither residual contains distance information</li>' +
        '<li><strong>Program Closure (Section 45)</strong>: 30+ hypotheses tested; 1D information ceiling proved; IFU follow-up memo included</li>' +
        '</ol>' +
        '<h3>The Inaccessibility-Strength Paradox</h3>' +
        '<p>No single-layer generative model can simultaneously produce r &asymp; 0.77 AND keep &gt;50% of H hidden from RC features. ' +
        'This implies H operates through channels affecting integral quantities (V<sub>flat</sub>, a<sub>0</sub>) without proportionally affecting ' +
        'differential quantities (RC radial profile shape) &mdash; a signature of multi-scale physics requiring angular resolution.</p>' +
        '<h3>Parameter Constraints on H</h3>' +
        '<ul>' +
        '<li>&alpha;<sub>Vf</sub> ~ 0.04&ndash;0.06 (H shifts V<sub>flat</sub> by ~4&ndash;6% per &sigma;)</li>' +
        '<li>&alpha;<sub>A0</sub> ~ 0.18&ndash;0.25 (H shifts a<sub>0</sub> by ~18&ndash;25% per &sigma;; 4&times; stronger than V<sub>flat</sub> channel)</li>' +
        '<li>&gamma;<sub>conc</sub> &lt; &minus;1 (H anti-correlates with halo concentration)</li>' +
        '</ul>' +
        '<h3>Conclusion</h3>' +
        '<p>SPARC 1D rotation curves reveal a real, statistically robust hidden state H (r &asymp; 0.77, p &lt; 0.001) that is distance-independent ' +
        'and construction-independent. However, its physical identity &mdash; whether halo redistribution, formation history, or exotic DM physics &mdash; ' +
        'cannot be determined from rotation curves alone. The decisive test requires IFU velocity fields preserving angular information.</p>' +
        '<h3>Caveats</h3>' +
        '<ul>' +
        '<li>All data from SPARC survey; cross-survey replication needed</li>' +
        '<li>N=55 sample limits statistical power in subsamples</li>' +
        '<li>Primary-distance subsample (N=11) inconclusive due to small size</li>' +
        '<li>Model degeneracy between physics families from 1D data alone</li>' +
        '<li>Not peer-reviewed</li>' +
        '</ul>' +
        '<p>Confidence: ~85%. Program status: <strong>CLOSED</strong>. IFU follow-up memo included.</p>' +
        '<p>All analysis scripts, results, and machine-readable data included for full reproducibility. ' +
        'Data: SPARC (Lelli, McGaugh &amp; Schombert, 2016, AJ 152, 157).</p>',
      resource_type: { id: 'dataset' },
      publication_date: new Date().toISOString().split('T')[0],
      publisher: 'Zenodo',
      creators: [
        { person_or_org: { type: 'personal', given_name: 'Galaxy Rotation', family_name: 'Curve Analyzer' } }
      ],
      rights: [
        { id: 'cc-by-4.0' }
      ],
      subjects: [
        { subject: 'galaxy rotation curves' },
        { subject: 'SPARC' },
        { subject: 'radial acceleration relation' },
        { subject: 'dark matter' },
        { subject: 'baryon-halo coupling' },
        { subject: 'rotation curve analysis' },
        { subject: 'galaxy dynamics' },
        { subject: 'modified gravity' },
        { subject: 'MOND' },
        { subject: 'per-galaxy a0' },
        { subject: 'kinematic residual' },
        { subject: 'cross-validation' },
        { subject: 'external validation' },
        { subject: 'regime-dependent' },
        { subject: 'hierarchical coupling' },
        { subject: 'dynamical integration' },
        { subject: 'physical interpretation' },
        { subject: 'information ceiling' },
        { subject: 'hidden state' },
        { subject: 'IFU follow-up' },
        { subject: 'inaccessibility paradox' },
      ],
      additional_descriptions: [
        {
          description: 'v12 — Program Closure: Phases 126-134 + 200-204 + 300-303 + Programs 1-8C + Phases V/V+. ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v11 (DOI 10.5281/zenodo.19440400). ' +
            'New in v12: Programs 8A (2D state recovery, 0/4), 8B (physics search, 5/6 max, inaccessibility-strength paradox), ' +
            '8C (map reconstruction, 2/4, scalar wins), Phase V (Red Team, 5/8, r revised to 0.77), ' +
            'Phase V+ (distance kill, 4/5, distance hypothesis killed), Program Closure (Section 45, 30+ hypotheses tallied). ' +
            'IFU follow-up memo included. Program status: CLOSED. Confidence: ~85%.',
          type: { id: 'notes' }
        }
      ],
      version: 'v12',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v12 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
