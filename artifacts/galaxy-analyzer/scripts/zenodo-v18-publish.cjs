#!/usr/bin/env node
const fs = require('fs');
const https = require('https');
const pathMod = require('path');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const CONCEPT_ID = '19446000';
const LATEST_PUBLISHED_ID = '19453653';
const ARCHIVE = pathMod.join(__dirname, '..', 'zenodo-archives', 'galaxy-rotation-curve-v18.tar.gz');
const BASE = 'zenodo.org';
const UA = 'GalaxyRotationCurveAnalyzer/1.0 (mailto:antminers99@gmail.com)';

function req(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE, path: urlPath, method,
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'Accept': 'application/json',
        'User-Agent': UA,
      },
    };
    if (body && contentType) {
      opts.headers['Content-Type'] = contentType;
      if (typeof body === 'string') opts.headers['Content-Length'] = Buffer.byteLength(body);
    }
    const r = https.request(opts, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('  HTTP ' + res.statusCode + ': ' + d.substring(0, 500));
          reject(new Error('HTTP ' + res.statusCode));
        } else {
          try { resolve(JSON.parse(d)); } catch { resolve(d); }
        }
      });
    });
    r.on('error', reject);
    if (body && typeof body !== 'string') body.pipe(r);
    else { if (body) r.write(body); r.end(); }
  });
}

function uploadToBucket(bucketUrl, localPath, remoteName) {
  return new Promise((resolve, reject) => {
    const stat = fs.statSync(localPath);
    const url = new URL(bucketUrl + '/' + remoteName);
    const opts = {
      hostname: url.hostname,
      path: url.pathname,
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'Content-Type': 'application/octet-stream',
        'Content-Length': stat.size,
        'User-Agent': UA,
      },
    };
    const r = https.request(opts, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('  Upload HTTP ' + res.statusCode + ': ' + d.substring(0, 300));
          reject(new Error('Upload fail'));
        } else {
          console.log('  Uploaded: ' + remoteName + ' (' + (stat.size / 1024 / 1024).toFixed(1) + ' MB)');
          resolve();
        }
      });
    });
    r.on('error', reject);
    fs.createReadStream(localPath).pipe(r);
  });
}

async function main() {
  console.log('=== Zenodo v18 Publish ===');
  console.log('DM-5C: Minimal carrier confirmed, phase frozen\n');

  if (!fs.existsSync(ARCHIVE)) {
    console.error('Archive not found: ' + ARCHIVE);
    process.exit(1);
  }
  console.log('Archive: ' + ARCHIVE);
  console.log('Size: ' + (fs.statSync(ARCHIVE).size / 1024 / 1024).toFixed(1) + ' MB\n');

  console.log('Step 1: Create new version from latest published...');
  const newVer = await req('POST', '/api/deposit/depositions/' + LATEST_PUBLISHED_ID + '/actions/newversion', '', 'application/json');
  const draftUrl = newVer.links?.latest_draft;
  if (!draftUrl) {
    console.error('No draft URL returned. Response: ' + JSON.stringify(newVer).substring(0, 300));
    process.exit(1);
  }
  const draftId = draftUrl.split('/').pop();
  console.log('  Draft created: ID=' + draftId);

  console.log('\nStep 2: Fetch draft details...');
  const draft = await req('GET', '/api/deposit/depositions/' + draftId, null, null);
  console.log('  State: ' + draft.state);
  const bucketUrl = draft.links?.bucket;
  console.log('  Bucket: ' + bucketUrl);

  console.log('\nStep 3: Delete old files...');
  const files = await req('GET', '/api/deposit/depositions/' + draftId + '/files', null, null);
  console.log('  Found ' + files.length + ' existing files');
  for (const f of files) {
    console.log('    Deleting: ' + f.filename);
    await req('DELETE', '/api/deposit/depositions/' + draftId + '/files/' + f.id, null, null);
  }

  console.log('\nStep 4: Upload new archive...');
  await uploadToBucket(bucketUrl, ARCHIVE, 'galaxy-rotation-curve-v18.tar.gz');

  console.log('\nStep 5: Update metadata...');
  const meta = {
    metadata: {
      title: 'Per-Galaxy a0 Variation: Angular Velocity-Field Complexity as Carrier of Hidden State (v18 — Programs 1-12, DM-5C Minimal Carrier Confirmed)',
      upload_type: 'dataset',
      publication_date: '2026-04-07',
      description:
        '<p><b>v18 — DM-5C complete: minimal carrier confirmed, phase frozen.</b></p>' +
        '<p>Single archive (1.3 MB): analysis scripts, result JSONs, source code, documentation.</p>' +
        '<p><b>Core Result (unchanged):</b> Hidden state H drives bilateral Vf_Resid-a0_Resid coupling at r = 0.77 (LOO, p &lt; 0.001, N = 55), replicated on N = 59 independent galaxies. 1D information ceiling: 88.5% inaccessible. Red team: 11/11 pass.</p>' +
        '<p><b>DM-5C: Extra Physics or Overfitting?</b></p>' +
        '<ul>' +
        '<li>M3 (3 params) collapses: LOO R&sup2; = &minus;1.055. The R&sup2; = 0.90 was overfitting.</li>' +
        '<li>shapeAmplitude is the minimal robust carrier of H (r = 0.80, LOO 7/7, bootstrap 97.8%).</li>' +
        '<li>outerSupport has independent information (partial r = 0.50), but N = 7 is insufficient to confirm.</li>' +
        '<li>downstreamQuietness is downstream, not causal (partial r = &minus;0.04 controlling for shape).</li>' +
        '</ul>' +
        '<p><b>Complete DM-5 Series:</b> DM-5A (CDM triaxial: 3/5 pass), DM-5B (hiddenness paradox RESOLVED), DM-5C (minimal carrier confirmed).</p>' +
        '<p><b>Model Status:</b> DEAD: CDM smooth, MOND, Fuzzy DM. No signal: SIDM, WDM. LEADING + SUFFICIENT: CDM + complex halo shape.</p>' +
        '<p><b>Phase Status:</b> FROZEN. No new analysis until new 2D kinematic data (WALLABY, MeerKAT).</p>' +
        '<p>Data: SPARC (Lelli, McGaugh &amp; Schombert 2016); THINGS (Walter et al. 2008).</p>',
      creators: [{ name: 'Fnd89', affiliation: 'Independent Researcher' }],
      access_right: 'open',
      license: 'cc-by-4.0',
      version: 'v18',
      keywords: ['galaxy rotation curves', 'dark matter', 'SPARC', 'THINGS', 'radial acceleration relation', 'halo shape', 'overfitting'],
    },
  };
  await req('PUT', '/api/deposit/depositions/' + draftId, JSON.stringify(meta), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 6: Publish...');
  const pub = await req('POST', '/api/deposit/depositions/' + draftId + '/actions/publish', '', 'application/json');
  console.log('\n=== PUBLISHED ===');
  console.log('DOI: ' + (pub.doi || 'check zenodo'));
  console.log('URL: ' + (pub.links?.html || pub.links?.record_html || 'check zenodo'));
  console.log('Concept DOI: 10.5281/zenodo.' + CONCEPT_ID);
}

main().catch(e => { console.error('\nFAILED: ' + e.message); process.exit(1); });
