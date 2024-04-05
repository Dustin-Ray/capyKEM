## capyKEM - research into quantum-resistant algorithm design

This repo is a pure rust, no-std interpretation of [FIPS 203 (draft)](https://csrc.nist.gov/pubs/fips/203/ipd) which leverages a module learning-with-errors based construction aimed towards providing a secure means of key exchange when faced with a potential quantum adversary.

### Current working items:

- [ ] document in style of FIPS
- [ ] support other two parameter sets
- [ ] parameterize sample_poly_cbd over eta
- [x] condense encoding/decoding to single function
- [x] parameterize encoding/decoding over d
- [ ] implement API-level functions
- [ ] replace usage of ```Vec``` with ```hybridarray```
