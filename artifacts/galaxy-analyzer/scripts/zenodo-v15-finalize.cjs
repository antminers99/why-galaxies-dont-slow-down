#!/usr/bin/env node
const fs = require('fs');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const DRAFT_ID = '19446498';
const ARCHIVE = '/tmp/zenodo-v15-package/galaxy-rotation-curve-v15.tar.gz';
const BASE = 'zenodo.org';

function req(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE, path: urlPath, method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'Accept': 'application/json' },
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
          console.error('  HTTP ' + res.statusCode + ': ' + d.substring(0, 300));
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

async function deleteFile(key) {
  await req('DELETE', '/api/records/' + DRAFT_ID + '/draft/files/' + encodeURIComponent(key), null, null);
}

async function uploadFile(localPath, remoteName) {
  const stat = fs.statSync(localPath);
  await req('POST', '/api/records/' + DRAFT_ID + '/draft/files',
    JSON.stringify([{ key: remoteName }]), 'application/json');
  await new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: '/api/records/' + DRAFT_ID + '/draft/files/' + encodeURIComponent(remoteName) + '/content',
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'Content-Type': 'application/octet-stream',
        'Content-Length': stat.size,
      },
    };
    const r = https.request(opts, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => {
        if (res.statusCode >= 400) { console.error('Upload err:', res.statusCode); reject(new Error('fail')); }
        else resolve();
      });
    });
    r.on('error', reject);
    fs.createReadStream(localPath).pipe(r);
  });
  await req('POST', '/api/records/' + DRAFT_ID + '/draft/files/' + encodeURIComponent(remoteName) + '/commit',
    '', 'application/json');
  console.log('  Uploaded: ' + remoteName + ' (' + stat.size + ' bytes)');
}

async function main() {
  console.log('=== Zenodo v15 Finalize: Single Archive ===\n');

  console.log('Step 1: Fetch draft...');
  const draft = await req('GET', '/api/records/' + DRAFT_ID + '/draft', null, null);
  console.log('  Title: ' + (draft.metadata?.title || 'unknown'));

  console.log('\nStep 2: Delete all existing files...');
  const files = await req('GET', '/api/records/' + DRAFT_ID + '/draft/files', null, null);
  const entries = files.entries || [];
  console.log('  Found ' + entries.length + ' files to delete');
  for (const f of entries) {
    await deleteFile(f.key);
  }
  console.log('  Done.');

  console.log('\nStep 3: Upload archive...');
  await uploadFile(ARCHIVE, 'galaxy-rotation-curve-v15.tar.gz');

  console.log('\nStep 4: Update metadata...');
  const meta = draft.metadata;
  meta.title = 'Per-Galaxy a0 Variation: Angular Velocity-Field Complexity as Carrier of Hidden State (v15 — Programs 1-11, post-submission verification)';
  meta.version = 'v15';
  meta.description =
    '<p><b>v15 — Post-submission interpretive refinement (Programs 1-11).</b></p>' +
    '<p>Single archive: 206 analysis scripts, 205 result JSONs, ' +
    '4 data files, 6 documents, CHANGELOG-v15.md.</p>' +
    '<p><b>Core Result (unchanged):</b> Hidden state H drives bilateral Vf_Resid-a0_Resid coupling ' +
    'at r = 0.77 (LOO, p &lt; 0.001, N = 55), replicated on N = 59 independent galaxies. ' +
    '1D information ceiling: 88.5% inaccessible. Red team: 11/11 pass (8/8 critical).</p>' +
    '<p><b>What changed:</b> Program 11 (16-bin azimuthal decomposition) shows m=3 carries ' +
    '58% of non-rotational power vs m=2 at 11%. H = total gravitational complexity beyond axisymmetry, ' +
    'not specifically m=2 triaxiality. CDM consistent but not uniquely favored.</p>' +
    '<p>MNRAS submission (v14) unchanged. See CHANGELOG-v15.md inside archive. ' +
    'Data: SPARC (Lelli, McGaugh &amp; Schombert 2016); THINGS (Walter et al. 2008).</p>';

  await req('PUT', '/api/records/' + DRAFT_ID + '/draft', JSON.stringify({ metadata: meta }), 'application/json');
  console.log('  Done.');

  console.log('\nStep 5: Publish...');
  const pub = await req('POST', '/api/records/' + DRAFT_ID + '/draft/actions/publish', '', 'application/json');
  console.log('\n=== PUBLISHED ===');
  console.log('DOI: ' + (pub.doi || 'check zenodo'));
  console.log('URL: ' + (pub.links?.self_html || pub.links?.html || 'check zenodo'));
}

main().catch(e => { console.error('\nFAILED: ' + e.message); process.exit(1); });
