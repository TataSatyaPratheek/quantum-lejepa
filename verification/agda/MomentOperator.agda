module MomentOperator where

data Nat : Set where
  zero : Nat
  suc  : Nat -> Nat
{-# BUILTIN NATURAL Nat #-}

data _==_ {A : Set} (a : A) : A -> Set where
  refl : a == a
{-# BUILTIN EQUALITY _==_ #-}

-- M₂ spectral decomposition: two blocks (Sym² and Λ²)
record M2Decomp : Set where
  field
    symDim  : Nat
    antiDim : Nat

-- Total dimension = symDim + antiDim
_+N_ : Nat -> Nat -> Nat
zero +N b = b
suc a +N b = suc (a +N b)

infixl 20 _+N_

-- The decomposition covers the full doubled space
totalDim : M2Decomp -> Nat
totalDim m = M2Decomp.symDim m +N M2Decomp.antiDim m
