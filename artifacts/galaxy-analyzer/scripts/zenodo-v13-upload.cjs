#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19444129';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/13.0', 'Accept': 'application/json' },
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
        'User-Agent': 'GalaxyAnalyzer/13.0',
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
  console.log('=== Zenodo v13 Upload: Carrier Identification + Red Team Verification ===');
  console.log('=== Programs 9 + 9V — Halo Triaxiality (m=2 mode) at ~95% confidence ===\n');

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
    { local: 'scripts/program9-phase901-crossmatch.cjs', remote: 'scripts/program9-phase901-crossmatch.cjs' },
    { local: 'scripts/program9-phase902-map-state.cjs', remote: 'scripts/program9-phase902-map-state.cjs' },
    { local: 'scripts/program9-phase903-decisive-test.cjs', remote: 'scripts/program9-phase903-decisive-test.cjs' },
    { local: 'scripts/program9-phase904-cosmo-search.cjs', remote: 'scripts/program9-phase904-cosmo-search.cjs' },
    { local: 'scripts/program9-phase905-carrier-id.cjs', remote: 'scripts/program9-phase905-carrier-id.cjs' },
    { local: 'scripts/program9v-red-team.cjs', remote: 'scripts/program9v-red-team.cjs' },
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
    { local: 'public/program9-phase901.json', remote: 'results/program9-phase901.json' },
    { local: 'public/program9-phase902.json', remote: 'results/program9-phase902.json' },
    { local: 'public/program9-phase903.json', remote: 'results/program9-phase903.json' },
    { local: 'public/program9-phase904.json', remote: 'results/program9-phase904.json' },
    { local: 'public/program9-phase905.json', remote: 'results/program9-phase905.json' },
    { local: 'public/program9v-red-team.json', remote: 'results/program9v-red-team.json' },
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
      title: 'Per-Galaxy a0 Variation: Carrier Identification via 2D Velocity Fields (v13 — Programs 1-9V, 52 sections)',
      description:
        '<p>v13 — <strong>Carrier identification + red team verification</strong>: the hidden state H driving bilateral VfResid&ndash;a0Resid coupling in SPARC rotation curves ' +
        'is identified as a robust outer-halo m=2 velocity-field mode, most consistent with halo triaxiality / oval distortion.</p>' +
        '<h3>Central Result</h3>' +
        '<p>A hidden common-cause state H drives bilateral coupling at <strong>r &asymp; 0.77</strong> (LOO, N=45). ' +
        'Program 9 identifies the observational carrier using real THINGS velocity field maps: the m=2 azimuthal mode correlates with bilateral excess at ' +
        '<strong>r = 0.847</strong> (p = 0.005, N=7). Gold pair NGC&nbsp;2841 vs NGC&nbsp;5055 shows 10.5&times; m=2 power difference despite near-identical baryonic mass. ' +
        'Cosmological Monte Carlo confirms triaxiality is <strong>required</strong> (concentration-only fails at 0% match rate).</p>' +
        '<h3>What is new in v13</h3>' +
        '<ol>' +
        '<li><strong>Program 9, Phase 901</strong>: IFU cross-match &mdash; 14/55 SPARC galaxies have 2D survey coverage; 4 matched pairs identified</li>' +
        '<li><strong>Phase 902</strong>: THINGS velocity field analysis &mdash; r(DQ, m2) = 0.899 in N=6; gold pair m2 ratio = 9.2&times;</li>' +
        '<li><strong>Phase 903</strong>: Decisive test &mdash; r(DQ, m2) = 0.847 (N=7), p = 0.005, LOO 7/7 positive, partial r|Vf = 0.750</li>' +
        '<li><strong>Phase 904</strong>: Cosmological MC &mdash; 8 models tested; 5 with triaxiality reproduce r &gt; 0.7; concentration-only FAILS</li>' +
        '<li><strong>Phase 905</strong>: Carrier identification &mdash; H1 (triaxiality) scores 23/23 (100%); H4 (exotic DM) scores 11/23 (48%)</li>' +
        '<li><strong>Program 9V</strong>: Red team verification &mdash; <strong>11/11 tests PASS</strong> (8/8 critical): independent pipeline replication (r = 0.911), ' +
        '19 sensitivity variations (all positive), bar contamination exclusion (signal in outer halo, r = 0.869), LOO 7/7, confounder control (partial r = 0.660), ' +
        'Program 8B reconciliation (consistent)</li>' +
        '<li><strong>Section 52</strong>: Conclusions &mdash; refined two-level claim (observational carrier ~95%, physical root cause ~90%)</li>' +
        '</ol>' +
        '<h3>The Refined Claim</h3>' +
        '<p>The hidden state H is strongly linked to a robust outer-halo m=2 velocity-field mode that survives independent reconstruction, sensitivity tests, ' +
        'bar exclusion, and confounder controls. The observational carrier is identified at ~95% confidence. ' +
        'The deeper physical origin is most consistent with a triaxial or closely related halo-shape state.</p>' +
        '<h3>Causal Chain</h3>' +
        '<p>Formation history &rarr; Halo triaxiality (b/a ~ 0.65&ndash;0.75) &rarr; m=2 velocity field power (angular, 1D-invisible) + enclosed mass shift (integral, 1D-visible) ' +
        '&rarr; Vflat &amp; a0 shifts (4:1 asymmetry) &rarr; bilateral coupling (r = 0.77 LOO)</p>' +
        '<h3>Caveats</h3>' +
        '<ul>' +
        '<li>N=7 2D sample (need N&ge;20 for definitive)</li>' +
        '<li>THINGS-only (cross-survey replication with MaNGA/CALIFA/PHANGS needed)</li>' +
        '<li>H1/H3 degeneracy (triaxiality vs assembly history as root cause)</li>' +
        '<li>Not peer-reviewed</li>' +
        '</ul>' +
        '<p>Confidence: <strong>~95%</strong> (observational carrier), <strong>~90%</strong> (physical root cause). ' +
        'Path to 99%: Program 10 (N&ge;20 multi-survey 2D replication).</p>' +
        '<p>All analysis scripts, results, and machine-readable data included for full reproducibility. ' +
        'Data: SPARC (Lelli, McGaugh &amp; Schombert, 2016, AJ 152, 157); THINGS (Walter et al. 2008, AJ 136, 2563).</p>',
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
        { subject: 'galaxy dynamics' },
        { subject: 'MOND' },
        { subject: 'per-galaxy a0' },
        { subject: 'hidden state' },
        { subject: 'halo triaxiality' },
        { subject: 'velocity field' },
        { subject: 'IFU kinematics' },
        { subject: 'THINGS survey' },
        { subject: 'm=2 azimuthal mode' },
        { subject: 'oval distortion' },
        { subject: 'information ceiling' },
        { subject: 'inaccessibility paradox' },
        { subject: 'carrier identification' },
        { subject: 'red team verification' },
        { subject: 'cross-validation' },
      ],
      additional_descriptions: [
        {
          description: 'v13 — Carrier Identification + Red Team: Programs 1-9V, 52 manuscript sections. ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v12 (DOI 10.5281/zenodo.19444129). ' +
            'New in v13: Program 9 (Phases 901-905, carrier identification from 2D THINGS velocity fields), ' +
            'Program 9V (11/11 red team tests PASS), Section 52 (conclusions). ' +
            'Carrier: outer-halo m=2 mode (~95% confidence), physical origin: halo triaxiality (~90%). ' +
            '81 files: 37 scripts, 34 results, 3 docs, 5 data files, 2 docs.',
          type: { id: 'notes' }
        }
      ],
      version: 'v13',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v13 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
