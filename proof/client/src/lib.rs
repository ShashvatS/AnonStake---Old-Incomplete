extern crate libc;
extern crate pairing;
extern crate rand;
extern crate num_cpus;
extern crate futures;
extern crate futures_cpupool;
extern crate bit_vec;
extern crate crossbeam;
extern crate byteorder;

use std::ffi::CStr;
use libc::c_char;

use std::fs::File;
use futures::Future;
use std::slice;
use std::mem::size_of;

use std::sync::Arc;
use std::vec::Vec;
use std::ptr;

pub mod multicore;
pub mod multiexp;
pub mod domain;
pub mod groth16;

use std::ops::{Add, Sub};
use std::fmt;
use std::error::Error;
use std::io;
use std::marker::PhantomData;

// For randomness (during paramgen and proof generation)
use rand::{thread_rng};
use rand::Rng;


// For benchmarking
use std::time::{Duration, Instant};

// Bring in some tools for using pairing-friendly curves
use pairing::{
    Engine,
    PrimeField,
    Field,
    CurveProjective,
    CurveAffine
};

use domain::Scalar;

// We're going to use the BLS12-381 pairing-friendly elliptic curve.
use pairing::bls12_381::{
    Bls12
};

// We're going to use the Groth16 proving system.
use groth16::{
    Proof,
    generate_random_parameters,
    ParameterSource,
    Parameters
};

use multiexp::{
    FullDensity,
    multiexp
};

use multicore::Worker;

//begin copy paste from bellman's lib file
//=============================================================================


/// Computations are expressed in terms of arithmetic circuits, in particular
/// rank-1 quadratic constraint systems. The `Circuit` trait represents a
/// circuit that can be synthesized. The `synthesize` method is called during
/// CRS generation and during proving.
pub trait Circuit<E: Engine> {
    /// Synthesize the circuit into a rank-1 quadratic constraint system
    fn synthesize<CS: ConstraintSystem<E>>(
        self,
        cs: &mut CS
    ) -> Result<(), SynthesisError>;
}

/// Represents a variable in our constraint system.
#[derive(Copy, Clone, Debug)]
pub struct Variable(Index);

impl Variable {
    /// This constructs a variable with an arbitrary index.
    /// Circuit implementations are not recommended to use this.
    pub fn new_unchecked(idx: Index) -> Variable {
        Variable(idx)
    }

    /// This returns the index underlying the variable.
    /// Circuit implementations are not recommended to use this.
    pub fn get_unchecked(&self) -> Index {
        self.0
    }
}

/// Represents the index of either an input variable or
/// auxillary variable.
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum Index {
    Input(usize),
    Aux(usize)
}

/// This represents a linear combination of some variables, with coefficients
/// in the scalar field of a pairing-friendly elliptic curve group.
#[derive(Clone)]
pub struct LinearCombination<E: Engine>(Vec<(Variable, E::Fr)>);

impl<E: Engine> AsRef<[(Variable, E::Fr)]> for LinearCombination<E> {
    fn as_ref(&self) -> &[(Variable, E::Fr)] {
        &self.0
    }
}

impl<E: Engine> LinearCombination<E> {
    pub fn zero() -> LinearCombination<E> {
        LinearCombination(vec![])
    }
}

impl<E: Engine> Add<(E::Fr, Variable)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(mut self, (coeff, var): (E::Fr, Variable)) -> LinearCombination<E> {
        self.0.push((var, coeff));

        self
    }
}

impl<E: Engine> Sub<(E::Fr, Variable)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, (mut coeff, var): (E::Fr, Variable)) -> LinearCombination<E> {
        coeff.negate();

        self + (coeff, var)
    }
}

impl<E: Engine> Add<Variable> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, other: Variable) -> LinearCombination<E> {
        self + (E::Fr::one(), other)
    }
}

impl<E: Engine> Sub<Variable> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, other: Variable) -> LinearCombination<E> {
        self - (E::Fr::one(), other)
    }
}

impl<'a, E: Engine> Add<&'a LinearCombination<E>> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(mut self, other: &'a LinearCombination<E>) -> LinearCombination<E> {
        for s in &other.0 {
            self = self + (s.1, s.0);
        }

        self
    }
}

impl<'a, E: Engine> Sub<&'a LinearCombination<E>> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(mut self, other: &'a LinearCombination<E>) -> LinearCombination<E> {
        for s in &other.0 {
            self = self - (s.1, s.0);
        }

        self
    }
}

impl<'a, E: Engine> Add<(E::Fr, &'a LinearCombination<E>)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(mut self, (coeff, other): (E::Fr, &'a LinearCombination<E>)) -> LinearCombination<E> {
        for s in &other.0 {
            let mut tmp = s.1;
            tmp.mul_assign(&coeff);
            self = self + (tmp, s.0);
        }

        self
    }
}

impl<'a, E: Engine> Sub<(E::Fr, &'a LinearCombination<E>)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(mut self, (coeff, other): (E::Fr, &'a LinearCombination<E>)) -> LinearCombination<E> {
        for s in &other.0 {
            let mut tmp = s.1;
            tmp.mul_assign(&coeff);
            self = self - (tmp, s.0);
        }

        self
    }
}

/// This is an error that could occur during circuit synthesis contexts,
/// such as CRS generation, proving or verification.
#[derive(Debug)]
pub enum SynthesisError {
    /// During synthesis, we lacked knowledge of a variable assignment.
    AssignmentMissing,
    /// During synthesis, we divided by zero.
    DivisionByZero,
    /// During synthesis, we constructed an unsatisfiable constraint system.
    Unsatisfiable,
    /// During synthesis, our polynomials ended up being too high of degree
    PolynomialDegreeTooLarge,
    /// During proof generation, we encountered an identity in the CRS
    UnexpectedIdentity,
    /// During proof generation, we encountered an I/O error with the CRS
    IoError(io::Error),
    /// During verification, our verifying key was malformed.
    MalformedVerifyingKey,
    /// During CRS generation, we observed an unconstrained auxillary variable
    UnconstrainedVariable
}

impl From<io::Error> for SynthesisError {
    fn from(e: io::Error) -> SynthesisError {
        SynthesisError::IoError(e)
    }
}

impl Error for SynthesisError {
    fn description(&self) -> &str {
        match *self {
            SynthesisError::AssignmentMissing => "an assignment for a variable could not be computed",
            SynthesisError::DivisionByZero => "division by zero",
            SynthesisError::Unsatisfiable => "unsatisfiable constraint system",
            SynthesisError::PolynomialDegreeTooLarge => "polynomial degree is too large",
            SynthesisError::UnexpectedIdentity => "encountered an identity element in the CRS",
            SynthesisError::IoError(_) => "encountered an I/O error",
            SynthesisError::MalformedVerifyingKey => "malformed verifying key",
            SynthesisError::UnconstrainedVariable => "auxillary variable was unconstrained"
        }
    }
}

impl fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        if let &SynthesisError::IoError(ref e) = self {
            write!(f, "I/O error: ")?;
            e.fmt(f)
        } else {
            write!(f, "{}", self.description())
        }
    }
}

/// Represents a constraint system which can have new variables
/// allocated and constrains between them formed.
pub trait ConstraintSystem<E: Engine>: Sized {
    /// Represents the type of the "root" of this constraint system
    /// so that nested namespaces can minimize indirection.
    type Root: ConstraintSystem<E>;

    /// Return the "one" input variable
    fn one() -> Variable {
        Variable::new_unchecked(Index::Input(0))
    }

    /// Allocate a private variable in the constraint system. The provided function is used to
    /// determine the assignment of the variable. The given `annotation` function is invoked
    /// in testing contexts in order to derive a unique name for this variable in the current
    /// namespace.
    fn alloc<F, A, AR>(
        &mut self,
        annotation: A,
        f: F
    ) -> Result<Variable, SynthesisError>
        where F: FnOnce() -> Result<E::Fr, SynthesisError>, A: FnOnce() -> AR, AR: Into<String>;

    /// Allocate a public variable in the constraint system. The provided function is used to
    /// determine the assignment of the variable.
    fn alloc_input<F, A, AR>(
        &mut self,
        annotation: A,
        f: F
    ) -> Result<Variable, SynthesisError>
        where F: FnOnce() -> Result<E::Fr, SynthesisError>, A: FnOnce() -> AR, AR: Into<String>;

    /// Enforce that `A` * `B` = `C`. The `annotation` function is invoked in testing contexts
    /// in order to derive a unique name for the constraint in the current namespace.
    fn enforce<A, AR, LA, LB, LC>(
        &mut self,
        annotation: A,
        a: LA,
        b: LB,
        c: LC
    )
        where A: FnOnce() -> AR, AR: Into<String>,
              LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
              LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
              LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>;

    /// Create a new (sub)namespace and enter into it. Not intended
    /// for downstream use; use `namespace` instead.
    fn push_namespace<NR, N>(&mut self, name_fn: N)
        where NR: Into<String>, N: FnOnce() -> NR;

    /// Exit out of the existing namespace. Not intended for
    /// downstream use; use `namespace` instead.
    fn pop_namespace(&mut self);

    /// Gets the "root" constraint system, bypassing the namespacing.
    /// Not intended for downstream use; use `namespace` instead.
    fn get_root(&mut self) -> &mut Self::Root;

    /// Begin a namespace for this constraint system.
    fn namespace<'a, NR, N>(
        &'a mut self,
        name_fn: N
    ) -> Namespace<'a, E, Self::Root>
        where NR: Into<String>, N: FnOnce() -> NR
    {
        self.get_root().push_namespace(name_fn);

        Namespace(self.get_root(), PhantomData)
    }
}

/// This is a "namespaced" constraint system which borrows a constraint system (pushing
/// a namespace context) and, when dropped, pops out of the namespace context.
pub struct Namespace<'a, E: Engine, CS: ConstraintSystem<E> + 'a>(&'a mut CS, PhantomData<E>);

impl<'cs, E: Engine, CS: ConstraintSystem<E>> ConstraintSystem<E> for Namespace<'cs, E, CS> {
    type Root = CS::Root;

    fn one() -> Variable {
        CS::one()
    }

    fn alloc<F, A, AR>(
        &mut self,
        annotation: A,
        f: F
    ) -> Result<Variable, SynthesisError>
        where F: FnOnce() -> Result<E::Fr, SynthesisError>, A: FnOnce() -> AR, AR: Into<String>
    {
        self.0.alloc(annotation, f)
    }

    fn alloc_input<F, A, AR>(
        &mut self,
        annotation: A,
        f: F
    ) -> Result<Variable, SynthesisError>
        where F: FnOnce() -> Result<E::Fr, SynthesisError>, A: FnOnce() -> AR, AR: Into<String>
    {
        self.0.alloc_input(annotation, f)
    }

    fn enforce<A, AR, LA, LB, LC>(
        &mut self,
        annotation: A,
        a: LA,
        b: LB,
        c: LC
    )
        where A: FnOnce() -> AR, AR: Into<String>,
              LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
              LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
              LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>
    {
        self.0.enforce(annotation, a, b, c)
    }

    // Downstream users who use `namespace` will never interact with these
    // functions and they will never be invoked because the namespace is
    // never a root constraint system.

    fn push_namespace<NR, N>(&mut self, _: N)
        where NR: Into<String>, N: FnOnce() -> NR
    {
        panic!("only the root's push_namespace should be called");
    }

    fn pop_namespace(&mut self)
    {
        panic!("only the root's pop_namespace should be called");
    }

    fn get_root(&mut self) -> &mut Self::Root
    {
        self.0.get_root()
    }
}

impl<'a, E: Engine, CS: ConstraintSystem<E>> Drop for Namespace<'a, E, CS> {
    fn drop(&mut self) {
        self.get_root().pop_namespace()
    }
}

/// Convenience implementation of ConstraintSystem<E> for mutable references to
/// constraint systems.
impl<'cs, E: Engine, CS: ConstraintSystem<E>> ConstraintSystem<E> for &'cs mut CS {
    type Root = CS::Root;

    fn one() -> Variable {
        CS::one()
    }

    fn alloc<F, A, AR>(
        &mut self,
        annotation: A,
        f: F
    ) -> Result<Variable, SynthesisError>
        where F: FnOnce() -> Result<E::Fr, SynthesisError>, A: FnOnce() -> AR, AR: Into<String>
    {
        (**self).alloc(annotation, f)
    }

    fn alloc_input<F, A, AR>(
        &mut self,
        annotation: A,
        f: F
    ) -> Result<Variable, SynthesisError>
        where F: FnOnce() -> Result<E::Fr, SynthesisError>, A: FnOnce() -> AR, AR: Into<String>
    {
        (**self).alloc_input(annotation, f)
    }

    fn enforce<A, AR, LA, LB, LC>(
        &mut self,
        annotation: A,
        a: LA,
        b: LB,
        c: LC
    )
        where A: FnOnce() -> AR, AR: Into<String>,
              LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
              LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
              LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>
    {
        (**self).enforce(annotation, a, b, c)
    }

    fn push_namespace<NR, N>(&mut self, name_fn: N)
        where NR: Into<String>, N: FnOnce() -> NR
    {
        (**self).push_namespace(name_fn)
    }

    fn pop_namespace(&mut self)
    {
        (**self).pop_namespace()
    }

    fn get_root(&mut self) -> &mut Self::Root
    {
        (**self).get_root()
    }
}

//=============================================================================

struct EmptyCircuit<E: Engine> {
    size: u32,
    x: Option<E::Fr>
}

impl<E: Engine> Circuit<E> for EmptyCircuit<E> {
    fn synthesize<CS: ConstraintSystem<E>>(
        self,
        cs: &mut CS
    ) -> Result<(), SynthesisError>
    {
        let x_value = self.x;
        for i in 0..self.size {
            let tmp = cs.alloc(||i.to_string(),|| {
                x_value.ok_or(SynthesisError::AssignmentMissing)
            })?;

            cs.enforce(
                || "dummy constraint",
                |lc| lc + tmp,
                |lc| lc + CS::one(),
                |lc| lc + tmp
            );
        }

        Ok(())
    }
}

fn create_params(size: u32, file: File) -> u32 {
    let params = {
        let c = EmptyCircuit::<Bls12> {
            size: size,
            x: None
        };

        generate_random_parameters(c, &mut thread_rng()).unwrap()
    };  

    match params.write(file) {
        Ok(_) => return 0,
        Err(_) => return 1,
    };
}

static mut MAYBEPARAMS: Option<Parameters<Bls12>>  = None;

fn load_params(check: bool, file: File) {
    unsafe {
        match Parameters::<Bls12>::read(file, check) {
            Ok(s) => MAYBEPARAMS = Some(s),
            Err(_) => MAYBEPARAMS = None
        }
    }
}

fn provePartial<E: Engine, P: ParameterSource<E>>(
    mut params: P,
    fft_precomp: Arc<Vec<<<E>::Fr as PrimeField>::Repr>>
) -> Result<Proof<E>, SynthesisError> {

    let worker = Worker::new();
    let vk = params.get_vk(fft_precomp.len())?;

    let h = multiexp(&worker, params.get_h(fft_precomp.len())?, FullDensity, fft_precomp);

    let mut g_a = vk.delta_g1.mul(E::Fr::zero());
    g_a.add_assign(&h.wait()?);

    let g_b = vk.delta_g2.mul(E::Fr::zero());

    Ok(Proof {
        a: g_a.into_affine(),
        b: g_b.into_affine(),
        c: g_a.into_affine()
    })
}

fn proveFull<E: Engine, P: ParameterSource<E>>(
    mut params: P,
    h: Arc<Vec<<<E>::Fr as PrimeField>::Repr>>,
    input_assignment: Arc<Vec<<<E>::Fr as PrimeField>::Repr>>,
    aux_assignment: Arc<Vec<<<E>::Fr as PrimeField>::Repr>>,
    r: E::Fr,
    s: E::Fr
) -> Result<Proof<E>, SynthesisError> {
    let worker = Worker::new();
    let vk = params.get_vk(h.len())?;

    let h = multiexp(&worker, params.get_h(h.len())?, FullDensity, h);
    let l = multiexp(&worker, params.get_l(aux_assignment.len())?, FullDensity, aux_assignment.clone());
    
    let (a_inputs_source, a_aux_source) = params.get_a(input_assignment.len(), aux_assignment.len())?;
    let a_inputs = multiexp(&worker, a_inputs_source, FullDensity, input_assignment.clone());
    let a_aux = multiexp(&worker, a_aux_source, FullDensity, aux_assignment.clone());

    let (b_g1_inputs_source, b_g1_aux_source) = params.get_b_g1(input_assignment.len(), aux_assignment.len())?;
    let b_g1_inputs = multiexp(&worker, b_g1_inputs_source, FullDensity, input_assignment.clone());
    let b_g1_aux = multiexp(&worker, b_g1_aux_source, FullDensity, aux_assignment.clone());

    let (b_g2_inputs_source, b_g2_aux_source) = params.get_b_g2(input_assignment.len(), aux_assignment.len())?;
    
    let b_g2_inputs = multiexp(&worker, b_g2_inputs_source, FullDensity, input_assignment);
    let b_g2_aux = multiexp(&worker, b_g2_aux_source, FullDensity, aux_assignment);

    if vk.delta_g1.is_zero() || vk.delta_g2.is_zero() {
        // If this element is zero, someone is trying to perform a
        // subversion-CRS attack.
        return Err(SynthesisError::UnexpectedIdentity);
    }

    let mut g_a = vk.delta_g1.mul(r);
    g_a.add_assign_mixed(&vk.alpha_g1);
    let mut g_b = vk.delta_g2.mul(s);
    g_b.add_assign_mixed(&vk.beta_g2);
    let mut g_c;
    {
        let mut rs = r;
        rs.mul_assign(&s);

        g_c = vk.delta_g1.mul(rs);
        g_c.add_assign(&vk.alpha_g1.mul(s));
        g_c.add_assign(&vk.beta_g1.mul(r));
    }
    let mut a_answer = a_inputs.wait()?;
    a_answer.add_assign(&a_aux.wait()?);
    g_a.add_assign(&a_answer);
    a_answer.mul_assign(s);
    g_c.add_assign(&a_answer);

    let mut b1_answer = b_g1_inputs.wait()?;
    b1_answer.add_assign(&b_g1_aux.wait()?);
    let mut b2_answer = b_g2_inputs.wait()?;
    b2_answer.add_assign(&b_g2_aux.wait()?);

    g_b.add_assign(&b2_answer);
    b1_answer.mul_assign(r);
    g_c.add_assign(&b1_answer);
    g_c.add_assign(&h.wait()?);
    g_c.add_assign(&l.wait()?);
 
    Ok(Proof {
        a: g_a.into_affine(),
        b: g_b.into_affine(),
        c: g_c.into_affine()
    })

}

unsafe fn ptr_to_arc(values: *mut u64, size: u32) -> Arc<Vec<<<Bls12 as Engine>::Fr as PrimeField>::Repr>>{
    let ptr = values as *const Scalar<Bls12>;
    let tmp = slice::from_raw_parts(ptr, size as usize);
    let a = Arc::new(tmp.into_iter().map(|s| s.0.into_repr()).collect::<Vec<_>>());
    return a;
}

#[no_mangle]
pub extern "C" fn hello_world() {
    println!("Hello world from rust!");
    println!("sizeof scalar {}", size_of::<Scalar<Bls12>>());
    println!("sizeof fieldr {}", size_of::<<Bls12 as Engine>::Fr>());
    println!("size of proof: {}", size_of::<Proof<Bls12>>());
}

#[no_mangle]
pub unsafe extern "C" fn initParams(size: u32, file: *const c_char) -> u32 {
    if file.is_null() {
        return 1;
    }

    let raw = CStr::from_ptr(file);

    let file_as_str = match raw.to_str() {
        Ok(s) => s,
        Err(_) => return 1,
    };

    let buffer = match File::create(file_as_str) {
        Ok(s) => s,
        Err(_) => return 2,
    };

    let result = create_params(size, buffer);
    if result == 0 {
        println!("Parameters created!");
    }
    else {
        println!("Error while creating parameters!");
    }
    

    return result;
}

#[no_mangle]
pub unsafe extern "C" fn loadParams(check: bool, file: *const c_char) -> u32 {
    if file.is_null() {
        return 1;
    }

    let raw = CStr::from_ptr(file);

    let file_as_str = match raw.to_str() {
        Ok(s) => s,
        Err(_) => return 1,
    };

    let buffer = match File::open(file_as_str) {
        Ok(s) => s,
        Err(_) => return 2,
    };

    load_params(check, buffer);
    println!("Parameters loaded!");

    return 0;
}

#[no_mangle]
pub unsafe extern "C" fn createProofPartial(size: u32, values: *mut u64) -> (*mut Proof<Bls12>, u32) {
    //let a = Arc::<Vec::<Bls12::Fr::Repr>>::new();
    let a = ptr_to_arc(values, size);
    println!("size of arc: {}", size_of::<Arc<Vec<<<Bls12 as Engine>::Fr as PrimeField>::Repr>>>());

    if let Some(params) = &MAYBEPARAMS {
        let proof = provePartial(params, a);
        println!("size of proof: {}", size_of::<Proof<Bls12>>());

        if let Ok(s) = proof {
            let boxedproof = Box::new(s);
            return (Box::into_raw(boxedproof), 0);
            //return (&mut s, 0);
        }
        else {
            return (ptr::null_mut(), 1);
        }
    }
    else {
        return (ptr::null_mut(), 1);
    }
}

#[no_mangle]
pub unsafe extern "C" fn createProofFull(size: u32, input_size: u32, aux_size: u32, values: *mut u64, input: *mut u64, aux: *mut u64) -> (*mut Proof<Bls12>, u32) {
    let h = ptr_to_arc(values, size);
    let inputs = ptr_to_arc(input, input_size);
    let aux_inputs = ptr_to_arc(aux, aux_size);

    let rng = &mut thread_rng();

    let r: <Bls12 as Engine>::Fr = rng.gen();
    let s: <Bls12 as Engine>::Fr = rng.gen();
    
    if let Some(params) = &MAYBEPARAMS {
        let proof = proveFull(params, h, inputs, aux_inputs, r, s);

        if let Ok(s) = proof {
            let boxedproof = Box::new(s);
            return (Box::into_raw(boxedproof), 0);
            //return (&mut s, 0);
        }
        else {
            return (ptr::null_mut(), 1);
        }
    }
    else {
        return (ptr::null_mut(), 1);
    }
}

