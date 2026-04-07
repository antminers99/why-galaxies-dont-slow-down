#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19446000';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/15.0', 'Accept': 'application/json' },
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
        'User-Agent': 'GalaxyAnalyzer/15.0',
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
  console.log('=== Zenodo v15 Upload: Post-Submission Verification (Program 11) ===');
  console.log('=== Interpretive Refinement: m=2 -> angular velocity-field complexity ===\n');

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
    { local: 'scripts/program10-multi-survey.cjs', remote: 'scripts/program10-multi-survey.cjs' },
    { local: 'scripts/program11-phase1-dm-landscape.cjs', remote: 'scripts/program11-phase1-dm-landscape.cjs' },
    { local: 'scripts/program11-testA-radial-gradient.cjs', remote: 'scripts/program11-testA-radial-gradient.cjs' },
    { local: 'scripts/program11-testB-mspectrum.cjs', remote: 'scripts/program11-testB-mspectrum.cjs' },
    { local: 'scripts/program11-phase2-verdict.cjs', remote: 'scripts/program11-phase2-verdict.cjs' },
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
    { local: 'public/replication/dark-matter-program/results/phase11-1-landscape.json', remote: 'results/program11-phase1-landscape.json' },
    { local: 'public/replication/dark-matter-program/results/testA-radial-gradient.json', remote: 'results/program11-testA-radial-gradient.json' },
    { local: 'public/replication/dark-matter-program/results/testB-mspectrum.json', remote: 'results/program11-testB-mspectrum.json' },
    { local: 'public/replication/dark-matter-program/results/phase11-2-verdict.json', remote: 'results/program11-phase2-verdict.json' },
    { local: 'public/phase56-frozen-baselines.json', remote: 'data/phase56-frozen-baselines.json' },
    { local: 'public/stage-A-master-table.json', remote: 'data/stage-A-master-table.json' },
    { local: 'public/sparc-table.json', remote: 'data/sparc-table.json' },
    { local: 'public/phase58a2-tidal-expansion.json', remote: 'data/phase58a2-tidal-expansion.json' },
    { local: 'public/sparc-results.json', remote: 'data/sparc-results.json' },
    { local: 'src/data/sparc-datasets.ts', remote: 'data/sparc-datasets.ts' },
    { local: 'public/replication/REPRODUCIBLE_RESULT.md', remote: 'REPRODUCIBLE_RESULT.md' },
    { local: 'public/replication/MANUSCRIPT-TRACK2.md', remote: 'MANUSCRIPT-TRACK2.md' },
    { local: 'public/replication/IFU-FOLLOWUP-MEMO.md', remote: 'IFU-FOLLOWUP-MEMO.md' },
    { local: 'public/replication/dark-matter-program/README.md', remote: 'program11/README.md' },
    { local: 'public/replication/dark-matter-program/CHANGELOG-v15.md', remote: 'CHANGELOG-v15.md' },
  ];

  let uploaded = 0, skipped = 0;
  for (const f of filesToUpload) {
    const fullPath = path.join(baseDir, f.local);
    if (fs.existsSync(fullPath)) {
      await uploadFile(draftId, fullPath, f.remote);
      uploaded++;
    } else {
      console.log('  SKIPPED (not found): ' + f.local);
      skipped++;
    }
  }
  console.log('\n  Uploaded: ' + uploaded + ', Skipped: ' + skipped);

  console.log('\nStep 4: Update metadata...');
  const metadata = {
    metadata: {
      title: 'Per-Galaxy a0 Variation: Angular Velocity-Field Complexity as Carrier of Hidden State (v15 — Programs 1-11, post-submission verification)',
      description:
        '<p><strong>v15 &mdash; Post-submission interpretive refinement.</strong> ' +
        'This version preserves the core hidden-state result but updates the interpretation of the 2D carrier ' +
        'from a uniquely isolated m=2 mode to a broader angular/non-axisymmetric velocity-field complexity, ' +
        'after additional post-submission verification (Program 11).</p>' +

        '<h3>Core Result (unchanged)</h3>' +
        '<p>A hidden common-cause state H drives bilateral coupling between V<sub>f,Resid</sub> and a<sub>0,Resid</sub> at ' +
        '<strong>r &asymp; 0.77</strong> (LOO, p &lt; 0.001, N=55, replicated on N=59 independent galaxies). ' +
        'The 1D information ceiling is genuine: 88.5% of H is inaccessible from rotation-curve features alone. ' +
        'Red team verification: 11/11 tests pass (8/8 critical).</p>' +

        '<h3>What changed in v15 (Program 11)</h3>' +
        '<p>Program 11 performed two post-submission verification tests using 16-bin azimuthal decomposition ' +
        '(vs. 8 bins in Program 9):</p>' +
        '<ul>' +
        '<li><strong>Test A (radial m=2 gradient)</strong>: All 7 THINGS galaxies show outer-dominant m=2, but this is ' +
        'a <em>universal</em> pattern (r(DQ, ratio) = &minus;0.26, p = 0.49) &mdash; most likely a velocity-field scaling effect, ' +
        'not evidence for SIDM core suppression.</li>' +
        '<li><strong>Test B (full m-spectrum, m=1&hellip;6)</strong>: Among non-rotational modes, m=3 carries <strong>58%</strong> of power ' +
        'vs. m=2 at only 11%. Correlations with DQ: r(DQ, m=3) = 0.786, r(DQ, m=5) = 0.802, ' +
        'r(DQ, m=2) = 0.516. The original r(DQ, m2) = 0.847 was partly inflated by azimuthal aliasing ' +
        '(8 bins, Nyquist limit &asymp; m=3).</li>' +
        '</ul>' +

        '<h3>Revised Interpretation</h3>' +
        '<p>H tracks <strong>total gravitational complexity beyond axisymmetry</strong> &mdash; not specifically halo triaxiality (m=2). ' +
        'The non-circular motions are dominated by odd modes (m=1, 3, 5), characteristic of spiral arm streaming, warps, ' +
        'and accretion-driven asymmetries. CDM remains consistent but is not uniquely favored. ' +
        'SIDM and Fuzzy DM are not specifically supported by these tests.</p>' +

        '<h3>What remains unchanged</h3>' +
        '<ol>' +
        '<li>Existence of hidden state H (bilateral coupling, r = 0.77)</li>' +
        '<li>1D information ceiling (88.5% inaccessible)</li>' +
        '<li>2D angular carrier exists (velocity fields carry H information)</li>' +
        '<li>Red team verification (11/11 pass)</li>' +
        '<li>Inaccessibility&ndash;strength paradox</li>' +
        '</ol>' +

        '<h3>MNRAS Submission Status</h3>' +
        '<p>The submitted manuscript (v14) remains unchanged at ScholarOne. ' +
        'This Zenodo version serves as the corrected reference. ' +
        'If reviewers raise the interpretation question, the refined analysis is ready.</p>' +

        '<h3>Discrimination Table</h3>' +
        '<table>' +
        '<tr><th>Test</th><th>CDM</th><th>SIDM</th><th>Fuzzy DM</th><th>Observed</th></tr>' +
        '<tr><td>A (radial gradient)</td><td>Flat</td><td>Core suppressed</td><td>N/A</td><td>Universal outer-dominant, no DQ dependence</td></tr>' +
        '<tr><td>B (m-spectrum)</td><td>m=2 dominant</td><td>N/A</td><td>Spread to m=3,4</td><td>m=3 dominant (58% non-rot), odd modes prevail</td></tr>' +
        '</table>' +

        '<p>All analysis scripts (42 files), results (38 files), and changelog included for full reproducibility. ' +
        'Data: SPARC (Lelli, McGaugh &amp; Schombert 2016); THINGS (Walter et al. 2008).</p>',
      resource_type: { id: 'dataset' },
      publication_date: new Date().toISOString().split('T')[0],
      publisher: 'Zenodo',
      creators: [
        { person_or_org: { type: 'personal', given_name: 'Fnd89', family_name: 'Independent Researcher' } }
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
        { subject: 'velocity field complexity' },
        { subject: 'non-axisymmetric structure' },
        { subject: 'IFU kinematics' },
        { subject: 'THINGS survey' },
        { subject: 'azimuthal modes' },
        { subject: 'm-spectrum' },
        { subject: 'information ceiling' },
        { subject: 'inaccessibility paradox' },
        { subject: 'carrier identification' },
        { subject: 'red team verification' },
        { subject: 'cross-validation' },
        { subject: 'post-submission verification' },
      ],
      additional_descriptions: [
        {
          description: 'v15 — Post-Submission Verification: Programs 1-11, interpretive refinement. ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v14 (DOI 10.5281/zenodo.19446000). ' +
            'New in v15: Program 11 (Phase 11.1 DM landscape survey, Test A radial m=2 gradient, ' +
            'Test B full m-spectrum m=1..6, Phase 11.2 combined verdict). ' +
            'Key change: Carrier interpretation updated from "m=2 halo triaxiality" to ' +
            '"angular velocity-field complexity (non-axisymmetric structure)" based on 16-bin ' +
            'azimuthal decomposition showing m=3 dominates non-rotational power (58% vs 11% for m=2). ' +
            'Core result (H detection, bilateral coupling, 1D ceiling) UNCHANGED. ' +
            'See CHANGELOG-v15.md for detailed before/after comparison.',
          type: { id: 'notes' }
        }
      ],
      version: 'v15',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v15 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
