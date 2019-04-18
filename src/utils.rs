use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Zero, One};

pub fn inverse_mod(mut a: BigInt, m: BigInt) -> BigInt {
    if a < BigInt::zero() || m <= a {
        a = a.mod_floor(&m);
    }

    let (mut c, mut d) = (a, m.clone());
    let (mut uc, mut vc, mut ud, mut vd) = (BigInt::one(), BigInt::zero(), BigInt::zero(), BigInt::one());
    while c != BigInt::zero() {
        let (q, r) = d.div_mod_floor(&c);
        d = c;
        c = r;

        println!("{}, {}, {}", q, c, d);

        let (nuc, nvc, nud, nvd) = (
            &ud - &q * &uc,
            &vd - &q * &vc,
            uc,
            vc
        );

        uc = nuc;
        vc = nvc;
        ud = nud;
        vd = nvd;
        println!("{}, {}, {}, {}", uc, vc, ud, vd);
    }

    assert_eq!(d, BigInt::one());
    if ud > BigInt::zero() {
        ud
    } else {
        ud + m
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_mod_test() {
        println!("For 5, 13");
        assert_eq!(inverse_mod(BigInt::from(5u32), BigInt::from(13u32)), BigInt::from(8u32));
        println!("For -5, 13");
        assert_eq!(inverse_mod(BigInt::from(-5i32), BigInt::from(13u32)), BigInt::from(5u32));
    }
}
