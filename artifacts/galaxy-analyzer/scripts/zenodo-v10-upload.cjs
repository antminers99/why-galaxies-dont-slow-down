#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19433840';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/10.0', 'Accept': 'application/json' },
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
        'User-Agent': 'GalaxyAnalyzer/10.0',
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
  console.log('=== Zenodo v10 Upload: Comprehensive External Validation ===');
  console.log('=== Phases 200-204: Complete External Validation Program ===\n');

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
    { local: 'public/stage-A-master-table.json', remote: 'data/stage-A-master-table.json' },
    { local: 'public/sparc-table.json', remote: 'data/sparc-table.json' },
    { local: 'public/phase58a2-tidal-expansion.json', remote: 'data/phase58a2-tidal-expansion.json' },
    { local: 'public/sparc-results.json', remote: 'data/sparc-results.json' },
    { local: 'public/replication/REPRODUCIBLE_RESULT.md', remote: 'REPRODUCIBLE_RESULT.md' },
    { local: 'public/replication/MANUSCRIPT-TRACK2.md', remote: 'MANUSCRIPT-TRACK2.md' },
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
      title: 'Comprehensively Validated Regime-Dependent Coupling Law for Per-Galaxy a0 Variation in SPARC Galaxies (v10 — Phases 126-134, 200-204)',
      description:
        '<p>v10 — Phases 126-134 (internal analysis) + Phases 200-204 (comprehensive external validation on N=59 independent SPARC galaxies).</p>' +
        '<h3>Central Result</h3>' +
        '<p>VfResid (kinematic residual) dominance over structural/environmental channels is <strong>independent of environmental proxy quality</strong>. ' +
        'The complete hierarchical coupling law replicates externally with all 8 hierarchy checks passing.</p>' +
        '<h3>External Validation Program (Phases 200-204)</h3>' +
        '<ol>' +
        '<li><strong>Phase 200-201</strong>: Data assembly (N=59) and blind prediction with frozen N=45 coefficients</li>' +
        '<li><strong>Phase 202</strong>: Complete hierarchy replication — 8/8 checks pass, VfResid dominates in all regimes, P&gt;99.6%</li>' +
        '<li><strong>Phase 203</strong>: Environmental decomposition — improved logMhost (r: 0.232&rarr;0.542) proves Core is real but VfResid still dominates by 30-76pp</li>' +
        '<li><strong>Phase 204</strong>: Final synthesis — lhOuter genuine but reduced after Mhost correction (+18.1pp&rarr;+8.3pp very-high-V)</li>' +
        '</ol>' +
        '<h3>Key External Numbers</h3>' +
        '<table>' +
        '<tr><th>Regime</th><th>N</th><th>VfResid-only gap</th><th>VfResid margin over improved Core</th></tr>' +
        '<tr><td>Full sample</td><td>59</td><td>+47.5%</td><td>+40.3 pp</td></tr>' +
        '<tr><td>Vflat &ge; 120</td><td>16</td><td>+48.5%</td><td>+30.0 pp</td></tr>' +
        '<tr><td>Vflat &ge; 180</td><td>8</td><td>+75.4%</td><td>+76.3 pp</td></tr>' +
        '<tr><td>Q=1 + Vflat &ge; 120</td><td>11</td><td>+63.8%</td><td>+43.0 pp</td></tr>' +
        '</table>' +
        '<h3>Key Numbers (Internal, N=45)</h3>' +
        '<ul>' +
        '<li>LOO gap: Core = 44.1%, Core+VfResid = 61.1%, 5-axis = 65.4%</li>' +
        '<li>Initial holdout (N=10): Core+VfResid = 56.9%, 5-axis = 66.0%, r = 0.801</li>' +
        '<li>VfResid drivers: haloK explains 29%, best-5 explains 52%, 36% irreducible</li>' +
        '</ul>' +
        '<h3>Caveats</h3>' +
        '<ul>' +
        '<li>All data from SPARC survey; cross-survey replication needed</li>' +
        '<li>High-Vflat external subsamples still small (N=8-16)</li>' +
        '<li>Best logMhost estimator (Model C) explains only ~29% of real variance</li>' +
        '<li>lhOuter contribution partially inflated by Mhost error absorption in initial estimates</li>' +
        '<li>Not peer-reviewed</li>' +
        '</ul>' +
        '<p>All analysis scripts, external validation program, and machine-readable results included for full reproducibility. ' +
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
        { subject: 'transfer learning' },
        { subject: 'environmental decomposition' },
      ],
      additional_descriptions: [
        {
          description: 'v10 — Phases 126-134 + 200-204 (Comprehensively Validated Coupling Law). ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v9 (DOI 10.5281/zenodo.19433840). ' +
            'New in v10: P202 Hierarchy Replication (8/8 checks pass, VfResid dominates all regimes, bootstrap P>99.6%), ' +
            'P203 logMhost Improvement (Model C r=0.542, Core gap -53% to -8.1%, VfResid still dominates by 30-76pp), ' +
            'P204 Final External Synthesis (lhOuter genuine but reduced: +18.1pp to +8.3pp very-high-V after Mhost correction). ' +
            'Central result: VfResid dominance is independent of environmental proxy quality. ' +
            'Updated MANUSCRIPT-TRACK2.md with Sections 3.8, revised 3.7/4/5/6.',
          type: { id: 'notes' }
        }
      ],
      version: 'v10',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v10 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
