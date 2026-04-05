#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19431868';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN },
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

function uploadFile(bucketUrl, filePath, fileName) {
  return new Promise((resolve, reject) => {
    const stat = fs.statSync(filePath);
    const parsedBucket = new URL(bucketUrl + '/' + fileName);
    const opts = {
      hostname: parsedBucket.hostname,
      path: parsedBucket.pathname,
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
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
          reject(new Error('Upload failed'));
        } else {
          console.log('  Uploaded: ' + fileName + ' (' + stat.size + ' bytes)');
          resolve();
        }
      });
    });
    req.on('error', reject);
    fs.createReadStream(filePath).pipe(req);
  });
}

async function main() {
  console.log('=== Zenodo v5 Upload: Phase 122 Writing Freeze ===\n');

  console.log('Step 1: Create new version from record ' + LATEST_ID + '...');
  const newVer = await zenodoRequest('POST', '/api/records/' + LATEST_ID + '/versions', '', 'application/json');
  const draftId = newVer.id;
  const bucketUrl = newVer.links.bucket;
  console.log('  Draft ID: ' + draftId);
  console.log('  Bucket: ' + bucketUrl);

  console.log('\nStep 2: Delete old files from draft...');
  const draftInfo = await zenodoRequest('GET', '/api/records/' + draftId + '/draft', null, null);
  if (draftInfo.files && draftInfo.files.length > 0) {
    for (const f of draftInfo.files) {
      try {
        await zenodoRequest('DELETE', '/api/records/' + draftId + '/draft/files/' + f.key, null, null);
        console.log('  Deleted: ' + f.key);
      } catch (e) {
        console.log('  Could not delete ' + f.key + ': ' + e.message);
      }
    }
  }

  console.log('\nStep 3: Upload new files...');
  const baseDir = path.join(__dirname, '..');
  const filesToUpload = [
    { local: 'public/replication/REPRODUCIBLE_RESULT.md', remote: 'REPRODUCIBLE_RESULT.md' },
    { local: 'public/replication/MANUSCRIPT.md', remote: 'MANUSCRIPT.md' },
    { local: 'public/replication/LITERATURE_SUMMARY.md', remote: 'LITERATURE_SUMMARY.md' },
    { local: 'public/replication/references.md', remote: 'references.md' },
    { local: 'scripts/phase120-mechanism-map.cjs', remote: 'scripts/phase120-mechanism-map.cjs' },
    { local: 'scripts/phase121-minimal-physical-model.cjs', remote: 'scripts/phase121-minimal-physical-model.cjs' },
    { local: 'scripts/phase117-anti-circularity.cjs', remote: 'scripts/phase117-anti-circularity.cjs' },
    { local: 'scripts/phase118-halo-anchoring.cjs', remote: 'scripts/phase118-halo-anchoring.cjs' },
    { local: 'scripts/phase119-model-competition.cjs', remote: 'scripts/phase119-model-competition.cjs' },
    { local: 'scripts/phase112-matched-falsification.cjs', remote: 'scripts/phase112-matched-falsification.cjs' },
    { local: 'scripts/phase110-uncertainty-propagation.cjs', remote: 'scripts/phase110-uncertainty-propagation.cjs' },
    { local: 'scripts/phase111-boundary-mapping.cjs', remote: 'scripts/phase111-boundary-mapping.cjs' },
    { local: 'scripts/phase108-mediation-dissection.cjs', remote: 'scripts/phase108-mediation-dissection.cjs' },
    { local: 'scripts/phase106-target-robustness.cjs', remote: 'scripts/phase106-target-robustness.cjs' },
    { local: 'scripts/phase105-death-match-sparse.cjs', remote: 'scripts/phase105-death-match-sparse.cjs' },
    { local: 'scripts/phase104-external-replication.cjs', remote: 'scripts/phase104-external-replication.cjs' },
    { local: 'scripts/phase102-residual-physics.cjs', remote: 'scripts/phase102-residual-physics.cjs' },
    { local: 'scripts/phase101-null-geometric-coupling.cjs', remote: 'scripts/phase101-null-geometric-coupling.cjs' },
    { local: 'public/phase120-mechanism-map.json', remote: 'results/phase120-mechanism-map.json' },
    { local: 'public/phase121-minimal-model.json', remote: 'results/phase121-minimal-model.json' },
    { local: 'public/phase117-anti-circularity.json', remote: 'results/phase117-anti-circularity.json' },
    { local: 'public/phase118-halo-anchoring.json', remote: 'results/phase118-halo-anchoring.json' },
    { local: 'public/phase119-model-competition.json', remote: 'results/phase119-model-competition.json' },
    { local: 'public/phase112-matched-falsification.json', remote: 'results/phase112-matched-falsification.json' },
  ];

  for (const f of filesToUpload) {
    const fullPath = path.join(baseDir, f.local);
    if (fs.existsSync(fullPath)) {
      await uploadFile(bucketUrl, fullPath, f.remote);
    } else {
      console.log('  SKIPPED (not found): ' + f.local);
    }
  }

  console.log('\nStep 4: Update metadata...');
  const metadata = {
    metadata: {
      title: 'Gas-to-Stellar State Is the Leading Predictor of Outer Support Requirement in High-Vflat SPARC Galaxies (v5 — Writing Freeze)',
      description: 'Phase 122 Writing Freeze. Analysis complete and frozen. ' +
        'Main finding: log(MHI/L3.6) is the strongest independent predictor of outer mass discrepancy ' +
        '(logOMD = log10 mean Vobs^2/Vbar^2 outer) in high-Vflat SPARC galaxies (Vflat >= 70 km/s, N=104). ' +
        'Adopted model: logOMD = 1.749 + 0.203*log(MHI/L3.6) - 0.101*log(Mbar), LOO R^2=0.584. ' +
        'Includes Phase 120 mechanism map (ratio beats both components), Phase 121 minimal physical model ' +
        '(2-var adopted, 3-var exploratory), anti-circularity battery (6/6 pass), halo anchoring (76% persists), ' +
        'model competition (no single model wins), and complete manuscript draft. ' +
        'All scripts and results included for full reproducibility. ' +
        'Data: SPARC (Lelli, McGaugh & Schombert, 2016, AJ 152, 157).',
      upload_type: 'dataset',
      publication_date: new Date().toISOString().split('T')[0],
      access_right: 'open',
      license: 'cc-by-4.0',
      keywords: [
        'galaxy rotation curves', 'SPARC', 'dark matter', 'gas fraction',
        'mass discrepancy', 'radial acceleration relation', 'gas-to-stellar balance',
        'outer support requirement', 'writing freeze',
      ],
      notes: 'v5 — Writing Freeze (Phase 122). Analysis FROZEN. ' +
        'No further phases unless external replication or falsification. ' +
        'Concept DOI: 10.5281/zenodo.19431363',
      version: 'v5',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v5 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
