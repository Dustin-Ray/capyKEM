use crate::{
    constants::parameter_sets::ParameterSet,
    math::{ntt::NttElement, ring_element::RingElement},
    Message,
};

impl<P: ParameterSet> Message<P> {
    fn k_pke_encrypt(&self, ek_pke: &[u8], r: &[u8; 32]) -> Vec<u8> {
        let k = P::K as usize;
        let mut n = 0;
        let t_hat = RingElement::from(&ek_pke[0..384 * k]);
        let rho: &[u8] = &ek_pke[384 * k..(384 * k) + 32];

        // Generate the matrix a_hat
        let mut a_hat = vec![NttElement::zero(); k * k];
        for i in 0..k {
            for j in 0..k {
                a_hat[i * k + j] = NttElement::sample_ntt(rho, i, j);
            }
        }

        // generate e1
        let mut e_1 = vec![RingElement::zero(); k];
        for e_elem in e_1.iter_mut().take(k) {
            *e_elem = RingElement::sample_poly_cbd::<P>(r, n);
            n += 1;
        }

        // generate e2
        let e2 = RingElement::sample_poly_cbd::<P>(r, n);

        // generate r, run ntt k times
        let mut r_hat = vec![NttElement::zero(); k];
        for r_elem in r_hat.iter_mut().take(k) {
            *r_elem = RingElement::sample_poly_cbd::<P>(r, n).into();
            n += 1;
        }

        // generate A transpose
        // TODO: this is likely not correct
        let mut u = vec![RingElement::zero(); k];
        for i in 0..u.len() {
            u[i] = e_1[i];
            for j in 0..r.len() {
                // Is addition happening in NTT domain as well?
                u[i] += (a_hat[i * k + j] * r_hat[j]).into() // into ring = NTT^(-1)
            }
        }

        // let mu = RingElement::byte_decode(&self.m).d;
        // TODO: NEED compress and decompress for ring

        todo!()
    }
}
