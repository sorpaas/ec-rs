mod utils;

use num_bigint::BigInt;
use num_traits::Zero;
use num_integer::Integer;

use core::ops::Add;
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

/// An elliptic curve, where y^2 = x^3 + a*x + b (mod p).
pub trait Curve: Clone + Eq + PartialEq + Debug {
    fn p() -> BigInt;
    fn a() -> BigInt;
    fn b() -> BigInt;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Eq, PartialEq, Debug)]
    /// Testing curve defined by y^2 = x^3 + 1x + 7
    struct TestCurve;

    impl Curve for TestCurve {
        fn p() -> BigInt { BigInt::from(13u32) }
        fn a() -> BigInt { BigInt::from(1u32) }
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
}
