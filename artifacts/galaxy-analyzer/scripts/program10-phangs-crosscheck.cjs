const fs = require('fs');
const path = require('path');
const https = require('https');

function httpGet(url) {
  return new Promise((resolve, reject) => {
    https.get(url, { timeout: 15000, headers: { 'User-Agent': 'GalaxyAnalyzer/1.0' } }, (res) => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        return httpGet(res.headers.location).then(resolve).catch(reject);
      }
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => resolve({ status: res.statusCode, data }));
    }).on('error', reject).on('timeout', function() { this.destroy(); reject(new Error('timeout')); });
  });
}

function httpGetBinary(url) {
  return new Promise((resolve, reject) => {
    https.get(url, { timeout: 60000, headers: { 'User-Agent': 'GalaxyAnalyzer/1.0' } }, (res) => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        return httpGetBinary(res.headers.location).then(resolve).catch(reject);
      }
      const chunks = [];
      res.on('data', chunk => chunks.push(chunk));
      res.on('end', () => resolve({ status: res.statusCode, data: Buffer.concat(chunks) }));
    }).on('error', reject).on('timeout', function() { this.destroy(); reject(new Error('timeout')); });
  });
}

const phangsTargets = [
  { name: 'NGC2903', things: true },
  { name: 'NGC3521', things: true },
];

async function tryDownloadPHANGS(name) {
  const lcName = name.toLowerCase().replace(/\s+/g, '');
  const baseUrl = 'https://www.canfar.net/storage/vault/list/phangs/RELEASES/PHANGS-ALMA/';

  console.log('\n  Attempting PHANGS-ALMA download for ' + name + '...');

  const momUrls = [
    'https://archive.stsci.edu/hlsps/phangs/alma/' + lcName + '/' + lcName + '_12m+7m+tp_co21_broad_mom1.fits',
    'https://www.canfar.net/storage/vault/files/phangs/RELEASES/PHANGS-ALMA/' + lcName + '/' + lcName + '_12m+7m+tp_co21_broad_mom1.fits',
  ];

  for (const url of momUrls) {
    console.log('    Trying: ' + url.split('/').slice(-2).join('/'));
    try {
      const resp = await httpGetBinary(url);
      if (resp.status === 200 && resp.data.length > 1000) {
        const outDir = path.join(__dirname, '..', 'data', 'multi-survey-2d', 'phangs');
        fs.mkdirSync(outDir, { recursive: true });
        const outFile = path.join(outDir, name + '_PHANGS_CO21_MOM1.fits');
        fs.writeFileSync(outFile, resp.data);
        console.log('    SUCCESS: ' + (resp.data.length / 1024).toFixed(0) + ' KB saved');
        return outFile;
      } else {
        console.log('    HTTP ' + resp.status + ' (' + resp.data.length + ' bytes)');
      }
    } catch(e) {
      console.log('    Error: ' + e.message);
    }
  }
  console.log('    NOT AVAILABLE via direct download');
  return null;
}

async function main() {
  console.log('='.repeat(72));
  console.log('PROGRAM 10.3 — PHANGS CO CROSS-CHECK');
  console.log('Independent tracer verification for THINGS galaxies');
  console.log('='.repeat(72));

  console.log('\n  Targets in both THINGS and PHANGS:');
  for (const t of phangsTargets) {
    console.log('    ' + t.name + (t.things ? ' (THINGS processed)' : ''));
  }

  const results = {};
  for (const t of phangsTargets) {
    const file = await tryDownloadPHANGS(t.name);
    results[t.name] = file ? 'downloaded' : 'not available';
  }

  console.log('\n\n  RESULTS:');
  for (const [name, status] of Object.entries(results)) {
    console.log('    ' + name + ': ' + status);
  }

  console.log('\n  NOTE: PHANGS CO data provides an INDEPENDENT TRACER test:');
  console.log('  If m=2 power shows same DQ correlation in CO as in HI,');
  console.log('  this rules out HI-specific systematic artifacts.');
  console.log('  This is a cross-check, not a sample expansion (same galaxies).');

  console.log('\n' + '='.repeat(72));
}

main().catch(e => console.error('Fatal:', e));
