require 'matrix'

class Matrix
  # Applies the matrix as a Moebius transformation to the complex
  # number z.
  def moebius(z)
    (self * Vector[z, 1]).inject(:/)
  end

  # Computes the matrix corresponding to the cross ratio of
  # three complex numbers
  def self.cross_ratio(zs)
    f = (zs[1] - zs[2]) / (zs[1] - zs[0])

    Matrix[[f, -f*zs[0]], [1, -zs[2]]]
  end

  alias_method :proj_inv, :adjugate
  alias_method :projective_inverse, :adjugate
end
