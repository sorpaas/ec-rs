mod utils;

use num_bigint::BigInt;
use num_traits::{Zero, One};
use num_integer::Integer;

use core::ops::{Add, Mul};
use core::marker::PhantomData;
use core::fmt::Debug;

/// A point value.
#[derive(Clone, Eq, PartialEq, Debug)]
pub enum PointValue {
    Infinity,
    Value {
        x: BigInt,
        y: BigInt,
    }
}

/// A point on a curve.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Point<C: Curve> {
    pub value: PointValue,
    _marker: PhantomData<C>,
}

impl<C: Curve> From<PointValue> for Point<C> {
    fn from(value: PointValue) -> Self {
        Self {
            value,
            _marker: PhantomData,
        }
    }
}

impl<C: Curve> Point<C> {
    pub fn value(self) -> Option<(BigInt, BigInt)> {
        match self.value {
            PointValue::Infinity => None,
            PointValue::Value { x, y } => Some((x, y)),
        }
    }

    pub fn infinity() -> Self {
        Self {
            value: PointValue::Infinity,
            _marker: PhantomData,
        }
    }

    pub fn is_infinity(&self) -> bool {
        match self.value {
            PointValue::Infinity => true,
            PointValue::Value { .. } => false,
        }
    }

    pub fn is_valid(&self) -> bool {
        match self.value {
            PointValue::Infinity => true,
            PointValue::Value { ref x, ref y } => {
                (y * y - (x * x * x + &C::a() * x + &C::b())).mod_floor(&C::p()) == BigInt::zero()
            },
        }
    }

    pub fn double(&self) -> Self {
        let (x, y) = match self.value {
            PointValue::Infinity => return Self::infinity(),
            PointValue::Value { ref x, ref y } => (x.clone(), y.clone()),
        };

        let l = ((BigInt::from(3u32) * &x * &x + C::a()) * utils::inverse_mod(BigInt::from(2u32) * &y, C::p())).mod_floor(&C::p());
        let x3 = (&l * &l - BigInt::from(2u32) * &x).mod_floor(&C::p());
        let y3 = (&l * (&x - &x3) - &y).mod_floor(&C::p());

        Self::from(PointValue::Value { x: x3, y: y3 })
    }
}

impl<C: Curve> Add for Point<C> {
    type Output = Point<C>;

    fn add(self, other: Point<C>) -> Point<C> {
        let (ox, oy) = match other.value {
            PointValue::Infinity => return self,
            PointValue::Value { ref x, ref y } => (x.clone(), y.clone()),
        };

        let (sx, sy) = match self.value {
            PointValue::Infinity => return other,
            PointValue::Value { ref x, ref y } => (x.clone(), y.clone()),
        };

        if sx == ox {
            return if (sy + oy).mod_floor(&C::p()) == BigInt::zero() {
                Point::infinity()
            } else {
                self.double()
            }
        }

        let l = ((&oy - &sy) * utils::inverse_mod(&ox - &sx, C::p())).mod_floor(&C::p());
        let x3 = (&l * &l - &sx - &ox).mod_floor(&C::p());
        let y3 = (&l * (&sx - &x3) - &sy).mod_floor(&C::p());

        Self::from(PointValue::Value { x: x3, y: y3 })
    }
}

impl<C: Curve> Mul<BigInt> for Point<C> {
    type Output = Point<C>;

    fn mul(self, mut other: BigInt) -> Point<C> {
        assert!(other >= BigInt::zero());

        if other == BigInt::zero() {
            Self::infinity()
        } else {
            let mut ret = self.clone();
            other -= BigInt::one();
            while other > BigInt::zero() {
                ret = ret + self.clone();
                other -= BigInt::one();
            }
            ret
        }
    }
}

/// An elliptic curve, where y^2 = x^3 + a*x + b (mod p).
pub trait Curve: Clone + Eq + PartialEq + Debug {
    fn p() -> BigInt;
    fn a() -> BigInt;
    fn b() -> BigInt;
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Num;

    #[derive(Clone, Eq, PartialEq, Debug)]
    /// Testing curve defined by y^2 = x^3 + 1x + 7
    struct TestCurve;

    impl Curve for TestCurve {
        fn p() -> BigInt { BigInt::from(13u32) }
        fn a() -> BigInt { BigInt::from(1u32) }
        fn b() -> BigInt { BigInt::from(7u32) }
    }

    #[derive(Clone, Eq, PartialEq, Debug)]
    /// secp256k1
    struct P256K1Curve;

    impl Curve for P256K1Curve {
        fn p() -> BigInt { BigInt::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16).unwrap() }
        fn a() -> BigInt { BigInt::zero() }
        fn b() -> BigInt { BigInt::from(7u32) }
    }

    #[test]
    fn point_addition() {
        let p1 = Point::<TestCurve>::from(PointValue::Value { x: BigInt::from(9u32), y: BigInt::from(11u32) });
        let p2 = Point::<TestCurve>::from(PointValue::Value { x: BigInt::from(4u32), y: BigInt::from(6u32) });
        assert!(p1.is_valid());
        assert!(p2.is_valid());

        let p3 = p1 + p2;

        assert_eq!(p3.value(), Some((BigInt::from(1u32), BigInt::from(10u32))));
    }

    #[test]
    fn secp256k1() {
        let g = Point::<P256K1Curve>::from(PointValue::Value { x: BigInt::from_str_radix("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16).unwrap(), y: BigInt::from_str_radix("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 16).unwrap() });

        assert_eq!(g.clone() * BigInt::one(), g);
        assert_eq!(
            g.clone() * BigInt::from(2u32),
            Point::<P256K1Curve>::from(PointValue::Value {
                x: BigInt::from_str_radix("C6047F9441ED7D6D3045406E95C07CD85C778E4B8CEF3CA7ABAC09B95C709EE5", 16).unwrap(),
                y: BigInt::from_str_radix("1AE168FEA63DC339A3C58419466CEAEEF7F632653266D0E1236431A950CFE52A", 16).unwrap()
            })
        );
        assert_eq!(
            g.clone() * BigInt::from(3u32),
            Point::<P256K1Curve>::from(PointValue::Value {
                x: BigInt::from_str_radix("F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9", 16).unwrap(),
                y: BigInt::from_str_radix("388F7B0F632DE8140FE337E62A37F3566500A99934C2231B6CB9FD7584B8E672", 16).unwrap()
            })
        );
    }
}
