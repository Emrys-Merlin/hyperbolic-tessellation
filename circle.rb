# frozen_string_literal: true
require_relative 'matrix.rb'

require 'cmath'
include CMath

# require 'gnuplotrb'
# include GnuplotRB

class Circle
  attr_reader :r, :c, :coeff

  def initialize(opts = {})
    if !opts[:coeff].nil?
      @coeff = opts[:coeff]
      set_r_and_c
    else
      @c = opts[:c] unless opts[:c].nil?
      @c = Complex(opts[:x], opts[:y]) unless opts[:x].nil? || opts[:y].nil?
      @r = opts[:r]
      set_coeff
    end
  end

  def set_coeff
    @coeff = if @r.nil?
               # line
               [0.0, @c * 1.0i, (@c * 1.0i).conjugate, 0.0]
             else
               # circle
               [1.0, -@c.conjugate, -@c, @c.abs**2 - @r**2]
             end
  end

  def set_r_and_c
    if @coeff[0].zero?
      # line
      @r = nil
      @c = @coeff[1] * -1.0i / @coeff[1].abs
    else
      # circle
      @c = -@coeff[2] / @coeff[0]
      @r = Math.sqrt((@coeff[1].abs / @coeff[0])**2 - @coeff[3] / @coeff[0])
    end
  end

  def plot_circle(color = 'blue')
    df = generalized_circle_points
    std_circle_plot(df, color)
  end

  def std_circle_plot(df, color = 'blue')
    [df, with: 'lines', lt: "rgb '#{color}'"]
  end

  def generalized_circle_points(step = 0.01)
    # line
    if @r.nil?

      v = @c / @c.abs

      x = []
      y = []

      (-1..1).step(step).each do |t|
        x << t * v.real
        y << t * v.imag
      end

    # circle
    else
      x = []
      y = []

      ratio = @r / @c.abs
      if ratio.abs <= 1
        alpha = Math.acos(ratio)

        gamma = @c.arg + PI

        (-alpha..alpha).step(step).each do |t|
          p = @c + @r * exp((t + gamma) * 1i)

          x << p.real
          y << p.imag
        end
      end
    end

    [x, y]
  end

  def get_triple
    # line
    if @r.nil?
      v = @c / @c.abs
      p1 = [v, 0, -v]

    # circle
    else
      alpha = Math.acos(@r / @c.abs)
      gamma = @c.arg + PI
      p1 = [@c + @r * exp((alpha + gamma) * 1i)]
      p1 << @c + @r * exp(gamma * 1i)
      p1 << @c + @r * exp((-alpha + gamma) * 1i)
    end

    p1
  end

  # Tests whether point p is contained in the open disk.
  # TODO: Handle lines
  def in_disk?(p)
    return false if @r.nil?
    (p - @c).abs < @r
  end

  # Tests if the two disks are contained in one another.
  # In that case it checks, whether both or none of the disks
  # contain p. In that case it returns 1 if "self" is closer to
  # p, if circ is closer to p it returns 2. In all other cases it
  # returns 1
  # Warning: So far this method only works for circles, not lines
  # TODO: handle lines
  def triple_inclusion(p, circ, eps = 0.001)
    # TODO: handle lines
    return 0 if @r.nil? || circ.r.nil?

    dist = (@c - circ.c).abs
    r_max = [@r, circ.r].max
    r_min = [@r, circ.r].min

    if dist + r_min >= r_max + eps
      return 0
    elsif in_disk?(p) && circ.in_disk?(p)
      return @r < circ.r ? 1 : 2
    elsif !in_disk?(p) && !circ.in_disk?(p)
      return @r > circ.r ? 1 : 2
    else
      return 0
    end
  end

  # Translate a generalized circle by a moebius transformation g
  # Not sure if correct.
  # TODO: Check Complex Analysis lecture notes
  def translate_by(g)
    a = @coeff[0]
    b = @coeff[1]
    c = @coeff[2]
    d = @coeff[3]

    h = g.proj_inv

    at = (a * h[0, 0].abs**2 + 2 * (b * h[0, 0] * h[1, 0].conjugate).real +
          d * h[1, 0].abs**2).real
    bt = a * h[0, 0] * h[0, 1].conjugate + b * h[0, 0] * h[1, 1].conjugate +
         c * h[0, 1].conjugate * h[1, 0] + d * h[1, 0] * h[1, 1].conjugate
    ct = bt.conjugate
    dt = (a * h[0, 1].abs**2 + 2 * (b * h[0, 1] * h[1, 1].conjugate).real +
          d * h[1, 1].abs**2).real

    puts 'Warning! Nonzeror interesect' if at.zero? && (dt != 0)
    Circle.new(coeff: [at, bt, ct, dt])
  end

  # Check whether the generalized circle is a geodesic in the poincare disk
  # model. If it is a line, some coefficients need to be zero. In the case of a circle
  # we check the orthogonoality via the pythagorean theorem.
  def geodesic?(eps = 0.01)
    if @coeff[0].zero?
      @coeff[3].abs < eps
    else
      (@c.abs**2 - 1.0 - @r**2)**2 < eps**2
    end
  end

  # Compute the corresponding Moebius transformation which
  # maps self to circ.
  def get_transformation(circ)
    p1 = get_triple
    p2 = circ.get_triple.reverse

    m1 = Matrix.cross_ratio(p1)
    m2 = Matrix.cross_ratio(p2).proj_inv

    m2 * m1
  end

  def self.get_bisector(p, q)
    raise 'Cannot compute bisector for two identical points.' if p == q

    qa = q.abs
    pa = p.abs

    a = qa**2 - pa**2
    b = -(1 - pa**2) * q.conjugate + (1 - qa**2) * p.conjugate
    c = b.conjugate
    d = a

    Circle.new(coeff: [a, b, c, d])
  end
end
