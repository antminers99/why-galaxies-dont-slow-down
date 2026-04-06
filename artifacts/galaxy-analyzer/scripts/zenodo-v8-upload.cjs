#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19433077';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/8.0', 'Accept': 'application/json' },
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
        'User-Agent': 'GalaxyAnalyzer/8.0',
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
  console.log('=== Zenodo v8 Upload: Transferable Hierarchical Coupling Law ===');
  console.log('=== Phases 133A-134: Regime Law + Second Channel + External Validation ===\n');

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
    { local: 'public/stage-A-master-table.json', remote: 'data/stage-A-master-table.json' },
    { local: 'public/sparc-table.json', remote: 'data/sparc-table.json' },
    { local: 'public/phase58a2-tidal-expansion.json', remote: 'data/phase58a2-tidal-expansion.json' },
    { local: 'public/sparc-results.json', remote: 'data/sparc-results.json' },
    { local: 'public/replication/REPRODUCIBLE_RESULT.md', remote: 'REPRODUCIBLE_RESULT.md' },
    { local: 'public/replication/MANUSCRIPT.md', remote: 'MANUSCRIPT.md' },
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
      title: 'Transferable Hierarchical Coupling Law for Per-Galaxy a0 Variation in SPARC Galaxies (v8 — Phases 126-134)',
      description:
        '<p>v8 — Phases 126-134: Transferable hierarchical coupling law for per-galaxy a0 variation, with external validation.</p>' +
        '<h3>Final Claim</h3>' +
        '<p>Per-galaxy a0 is governed by a <strong>transferable, regime-dependent, hierarchical baryon-halo coupling law</strong>:</p>' +
        '<ul>' +
        '<li>A 3-axis structural core (MHI, Mhost, MeanRun) sets the baseline</li>' +
        '<li>A dominant kinematic coupling channel (VfResid = residual kinematic excess beyond baryonic prediction) mediates halo physics into a0</li>' +
        '<li>A secondary outer-halo channel (lh_outerImprove) contributes independently</li>' +
        '<li>The coupling activates sharply at Vflat ~ 181 km/s; strongest in high-Vflat galaxies</li>' +
        '<li>~36% of the coupling signal is irreducibly encoded in multi-scale dynamical integration</li>' +
        '</ul>' +
        '<h3>External Validation (Phase 134)</h3>' +
        '<p>Trained on N=45 published-quality galaxies, tested on N=10 crude-quality holdout galaxies:</p>' +
        '<ul>' +
        '<li>Core alone: <strong>fails</strong> (transfer gap = -11.5%)</li>' +
        '<li>Core + raw Vflat: modest (gap = 10.0%)</li>' +
        '<li>Core + VfResid: <strong>strong</strong> (gap = 56.9%, r = 0.801)</li>' +
        '<li>5-axis (Core + VfResid + lhOuter): <strong>strongest</strong> (gap = 66.0%)</li>' +
        '</ul>' +
        '<p>All coefficient signs preserved. Regime-restricted training outperforms full training. The transferable signal is the residual dynamical excess channel, not raw structural variables.</p>' +
        '<h3>Key Numbers</h3>' +
        '<ul>' +
        '<li>Internal LOO gap (N=45): Core = 44.1%, Core+VfResid = 61.1%, 5-axis = 65.4%</li>' +
        '<li>External transfer gap (N=10 holdout): Core+VfResid = 56.9%, 5-axis = 66.0%</li>' +
        '<li>VfResid drivers: haloK explains 29%, best-5 explains 52%, 36% irreducible</li>' +
        '<li>Regime: VfResid delta = +27.8pp in Vflat&ge;120 regime; unstable in low-Vflat</li>' +
        '</ul>' +
        '<p>Caveats: holdout sample N=10 (directional, not conclusive); all crude quality; all Vflat&ge;120 (regime split not independently tested externally).</p>' +
        '<p>All analysis scripts and machine-readable results included for full reproducibility. Data: SPARC (Lelli, McGaugh &amp; Schombert, 2016, AJ 152, 157).</p>',
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
        { subject: 'state law' },
        { subject: 'kinematic residual' },
        { subject: 'cross-validation' },
        { subject: 'model selection' },
        { subject: 'external validation' },
        { subject: 'transfer learning' },
        { subject: 'hierarchical coupling' },
      ],
      additional_descriptions: [
        {
          description: 'v8 — Phases 126-134 (Transferable Hierarchical Coupling Law + External Validation). ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v7 (DOI 10.5281/zenodo.19433077). ' +
            'New in v8: P132A Halo Death Match (VFRESID_DOMINANT), P132B Mediation/Causal (VFRESID_DOMINANT_CHANNEL), ' +
            'P132C External Robustness (MOSTLY_ROBUST), P133A Regime Law (VERY_STRONG_REGIME + SHARP_TRANSITION at Vflat~181), ' +
            'P133B Second Channel (GENUINE_5TH_AXIS = lh_outerImprove), P133C Coupling Drivers (PARTIAL_COUPLING_LAW + DEEP_RESIDUAL), ' +
            'P134 External Validation (STRONG_TRANSFER, r=0.801 on N=10 holdout).',
          type: { id: 'notes' }
        }
      ],
      version: 'v8',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v8 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
