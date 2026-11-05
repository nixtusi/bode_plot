#!/usr/bin/env node

// A small App Store Connect API client. Credentials stay in the ignored local
// configuration and the JWT is kept in memory only.
import { createSign } from 'node:crypto';
import { readFile, stat } from 'node:fs/promises';
import https from 'node:https';
import { fileURLToPath } from 'node:url';

const apiOrigin = 'https://api.appstoreconnect.apple.com';
const configPath = fileURLToPath(new URL('../.appstoreconnect/config.json', import.meta.url));
const maxBodyBytes = 10 * 1024 * 1024;

function usage() {
  return `Usage: node scripts/app_store_connect.mjs [--method GET|POST|PATCH|DELETE] [--scope-read] [--data-file FILE] PATH\n\n` +
    `PATH must be an App Store Connect API path or an https://api.appstoreconnect.apple.com URL.\n` +
    `GET is the default. POST and PATCH require a JSON file via --data-file.\n`;
}

function parseArguments(args) {
  let method = 'GET';
  let dataFile;
  let path;
  let methodSpecified = false;
  let scopeRead = false;
  for (let i = 0; i < args.length; i++) {
    const arg = args[i];
    if (arg === '--help' || arg === '-h') return { help: true };
    if (arg === '--scope-read') {
      if (scopeRead) throw new Error('Invalid --scope-read.');
      scopeRead = true;
    } else if (arg === '--method') {
      if (methodSpecified || !args[i + 1]) throw new Error('Invalid --method.');
      method = args[++i].toUpperCase();
      methodSpecified = true;
    } else if (arg === '--data-file') {
      if (dataFile !== undefined || !args[i + 1]) throw new Error('Invalid --data-file.');
      dataFile = args[++i];
    } else if (arg.startsWith('-') || path !== undefined) {
      throw new Error('Invalid arguments.');
    } else {
      path = arg;
    }
  }
  if (!['GET', 'POST', 'PATCH', 'DELETE'].includes(method) || !path) {
    throw new Error('A supported method and API path are required.');
  }
  if ((method === 'POST' || method === 'PATCH') !== (dataFile !== undefined)) {
    throw new Error('Use --data-file only with POST or PATCH, where it is required.');
  }
  if (scopeRead && method !== 'GET') throw new Error('--scope-read is only available for GET requests.');
  return { method, dataFile, path, scopeRead };
}

function apiURL(path) {
  if (/[\x00-\x20\x7f]/.test(path) || (path.startsWith('/') && path.startsWith('//'))) {
    throw new Error('Invalid API path.');
  }
  let url;
  try {
    url = new URL(path, apiOrigin);
  } catch {
    throw new Error('Invalid API path.');
  }
  if (url.origin !== apiOrigin || url.username || url.password || url.hash || !url.pathname.startsWith('/v')) {
    throw new Error('Only App Store Connect API URLs are allowed.');
  }
  return url;
}

async function loadCredentials() {
  let config;
  try {
    config = JSON.parse(await readFile(configPath, 'utf8'));
  } catch {
    throw new Error('Could not read the local App Store Connect configuration.');
  }
  if (!config || typeof config.keyId !== 'string' || !config.keyId ||
      typeof config.issuerId !== 'string' || !config.issuerId ||
      typeof config.privateKeyPath !== 'string' || !config.privateKeyPath) {
    throw new Error('Local configuration must contain keyId, issuerId, and privateKeyPath.');
  }
  let privateKey;
  try {
    privateKey = await readFile(config.privateKeyPath, 'utf8');
  } catch {
    throw new Error('Could not read the configured private key.');
  }
  return { ...config, privateKey };
}

function makeToken({ keyId, issuerId, privateKey }, url, scopeRead) {
  const encode = value => Buffer.from(JSON.stringify(value)).toString('base64url');
  const now = Math.floor(Date.now() / 1000);
  const payload = { iss: issuerId, iat: now, exp: now + 600, aud: 'appstoreconnect-v1' };
  if (scopeRead) payload.scope = [`GET ${url.pathname}${url.search}`];
  const contents = `${encode({ alg: 'ES256', kid: keyId, typ: 'JWT' })}.${encode(payload)}`;
  try {
    const signer = createSign('SHA256');
    signer.update(contents);
    signer.end();
    const signature = signer.sign({ key: privateKey, dsaEncoding: 'ieee-p1363' }).toString('base64url');
    return `${contents}.${signature}`;
  } catch {
    throw new Error('Could not sign the API request with the configured private key.');
  }
}

async function loadBody(dataFile) {
  if (dataFile === undefined) return undefined;
  let body;
  try {
    if ((await stat(dataFile)).size > maxBodyBytes) throw new Error();
    body = await readFile(dataFile);
    if (body.length > maxBodyBytes) throw new Error();
    JSON.parse(body.toString('utf8'));
  } catch {
    throw new Error('Could not read a valid JSON request body (maximum 10 MiB).');
  }
  return body;
}

function request(url, method, body, token) {
  return new Promise((resolve, reject) => {
    const headers = { Authorization: `Bearer ${token}`, Accept: 'application/json' };
    if (body) {
      headers['Content-Type'] = 'application/json';
      headers['Content-Length'] = body.length;
    }
    const req = https.request({ protocol: 'https:', hostname: 'api.appstoreconnect.apple.com',
      method, path: `${url.pathname}${url.search}`, headers, timeout: 30000 }, response => {
      const chunks = [];
      let size = 0;
      response.on('data', chunk => {
        size += chunk.length;
        if (size > maxBodyBytes) {
          response.destroy();
          reject(new Error('API response exceeded 10 MiB.'));
        } else chunks.push(chunk);
      });
      response.on('end', () => resolve({ status: response.statusCode, body: Buffer.concat(chunks) }));
      response.on('error', () => reject(new Error('Could not read the API response.')));
    });
    req.on('timeout', () => req.destroy());
    req.on('error', () => reject(new Error('App Store Connect API request failed.')));
    if (body) req.write(body);
    req.end();
  });
}

function redact(value, token, privateKey) {
  if (typeof value === 'string') return value.replaceAll(token, '[REDACTED]').replaceAll(privateKey, '[REDACTED]');
  if (Array.isArray(value)) return value.map(item => redact(item, token, privateKey));
  if (value && typeof value === 'object') {
    return Object.fromEntries(Object.entries(value).map(([key, item]) => [key,
      /private.?key|secret|access.?token|authorization|password/i.test(key) ? '[REDACTED]' : redact(item, token, privateKey)]));
  }
  return value;
}

async function main() {
  const args = parseArguments(process.argv.slice(2));
  if (args.help) { process.stdout.write(usage()); return; }
  const url = apiURL(args.path);
  const body = await loadBody(args.dataFile);
  const credentials = await loadCredentials();
  const token = makeToken(credentials, url, args.scopeRead);
  const response = await request(url, args.method, body, token);
  if (response.status >= 300 && response.status < 400) throw new Error(`API returned HTTP ${response.status}; redirects are not followed.`);
  if (response.status < 200 || response.status >= 300) {
    let detail = '';
    try {
      const error = JSON.parse(response.body.toString('utf8')).errors?.[0];
      if (typeof error?.code === 'string') detail += ` ${error.code}`;
      if (typeof error?.detail === 'string') detail += `: ${error.detail}`;
    } catch { /* The HTTP status is still useful if the response is not JSON. */ }
    throw new Error(`API returned HTTP ${response.status}.${redact(detail, token, credentials.privateKey)}`);
  }
  if (response.body.length) {
    let json;
    try { json = JSON.parse(response.body.toString('utf8')); }
    catch { throw new Error('API returned a non-JSON response.'); }
    process.stdout.write(`${JSON.stringify(redact(json, token, credentials.privateKey), null, 2)}\n`);
  }
}

main().catch(error => {
  process.stderr.write(`${error.message}\n`);
  process.exitCode = 1;
});
