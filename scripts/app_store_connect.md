# App Store Connect API

This client uses Node.js built-in modules and a locally configured **Team API Key**. Its configuration is `.appstoreconnect/config.json`, which contains `keyId`, `issuerId`, and `privateKeyPath`. Keep that directory ignored by Git, and keep the `.p8` private key outside the repository. Do not share the configuration, private key, or API output containing sensitive account data.

Read an endpoint (GET is the default). `--scope-read` limits the short-lived token to that one GET request:

```sh
node scripts/app_store_connect.mjs --scope-read '/v1/apps?limit=1'
```

The API response is printed as JSON. An absolute pagination URL can be passed as the path only when it uses `https://api.appstoreconnect.apple.com`. Redirects are not followed. The client redacts recognizable credential fields and the token from JSON output, but API responses can still contain private account data.

For an explicitly authorized change, supply a JSON body from a local file for POST or PATCH, or use DELETE without a body:

```sh
node scripts/app_store_connect.mjs --method PATCH --data-file /path/to/request.json '/v1/resource/ID'
node scripts/app_store_connect.mjs --method DELETE '/v1/resource/ID'
```

The API key's App Store Connect role determines which endpoints are available. The client does not initiate changes unless a write method is passed.

Apple's documentation: [Generating Tokens for API Requests](https://developer.apple.com/documentation/appstoreconnectapi/generating-tokens-for-api-requests) and [App Store Connect API](https://developer.apple.com/documentation/appstoreconnectapi).
