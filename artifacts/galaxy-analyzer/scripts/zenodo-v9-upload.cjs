#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19433329';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/9.0', 'Accept': 'application/json' },
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
        'User-Agent': 'GalaxyAnalyzer/9.0',
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
  console.log('=== Zenodo v9 Upload: External Validation of Coupling Law ===');
  console.log('=== Phases 200-201: Broader External Validation (N=59) ===\n');

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
      title: 'Externally Validated Regime-Dependent Coupling Law for Per-Galaxy a0 Variation in SPARC Galaxies (v9 — Phases 126-134, 200-201)',
      description:
        '<p>v9 — Phases 126-134 (internal analysis) + Phases 200-201 (external validation on N=59 SPARC galaxies outside training sample).</p>' +
        '<h3>Central Result</h3>' +
        '<p>The coupling-law picture now has <strong>external support beyond the original N=45 sample</strong>. ' +
        'The strongest transfer appears in higher-Vflat galaxies and collapses in lower-Vflat systems, ' +
        'consistent with a <strong>regime-dependent law rather than a universal one</strong>.</p>' +
        '<h3>Hierarchical Structure (reproduced externally)</h3>' +
        '<ol>' +
        '<li>Structural core alone <strong>fails</strong> in all regimes (gap = -53% full sample)</li>' +
        '<li>VfResid is the <strong>dominant transferable channel</strong> (carries all positive gap)</li>' +
        '<li>Signal strengthens monotonically with Vflat</li>' +
        '<li>Low-Vflat regime fails, exactly as predicted by internal regime law</li>' +
        '</ol>' +
        '<h3>External Validation Summary (Phase 201)</h3>' +
        '<table>' +
        '<tr><th>Regime</th><th>N</th><th>Core+VfResid gap</th><th>r(VfResid, a&sub0;)</th></tr>' +
        '<tr><td>Full sample</td><td>59</td><td>+8.2%</td><td>0.713</td></tr>' +
        '<tr><td>Vflat &ge; 120</td><td>16</td><td>+34.3%</td><td>0.841</td></tr>' +
        '<tr><td>Vflat &ge; 180</td><td>8</td><td>+48.7%</td><td>0.858</td></tr>' +
        '<tr><td>Q=1 + Vflat &ge; 120</td><td>11</td><td><strong>+59.0%</strong></td><td><strong>0.830</strong></td></tr>' +
        '<tr><td>Vflat &lt; 120</td><td>43</td><td>-3.9%</td><td>&mdash;</td></tr>' +
        '</table>' +
        '<h3>Key Numbers (Internal)</h3>' +
        '<ul>' +
        '<li>Internal LOO gap (N=45): Core = 44.1%, Core+VfResid = 61.1%, 5-axis = 65.4%</li>' +
        '<li>Initial holdout (N=10): Core+VfResid = 56.9%, 5-axis = 66.0%, r = 0.801</li>' +
        '<li>VfResid drivers: haloK explains 29%, best-5 explains 52%, 36% irreducible</li>' +
        '</ul>' +
        '<h3>Caveats</h3>' +
        '<ul>' +
        '<li>External sample uses cruder logA0 (fewer RAR points) and estimated logMhost</li>' +
        '<li>High-Vflat external subsample still small (N=8-16)</li>' +
        '<li>All data from SPARC survey; cross-survey replication needed</li>' +
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
      ],
      additional_descriptions: [
        {
          description: 'v9 — Phases 126-134 + 200-201 (Externally Validated Regime-Dependent Coupling Law). ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v8 (DOI 10.5281/zenodo.19433329). ' +
            'New in v9: P200 Data Assembly (N=59 external SPARC galaxies with per-galaxy logA0, VfResid, logMhost, lhOuter), ' +
            'P201 Blind Prediction (frozen N=45 coefficients tested on N=59 external sample: MODERATE_TRANSFER overall +8.2%, ' +
            'HIGH_VFLAT_STRONG +34.3% to +59.0%, LOW_VFLAT_FAILS -3.9%, hierarchy preserved). ' +
            'Updated MANUSCRIPT-TRACK2.md with external validation section and revised conclusions.',
          type: { id: 'notes' }
        }
      ],
      version: 'v9',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v9 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
