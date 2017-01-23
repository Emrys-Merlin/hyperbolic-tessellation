require_relative '../circle.rb'
require 'minitest/autorun'
require 'minitest/pride'

class TestCircle < MiniTest::Test
  def setup
    @c = Circle.new(x: Math::sqrt(2), y: 0.0, r: 1.0)
    @prng = Random.new
  end

  # test that if a gc is given by radius and center that the
  # coefficients are computed correctly
  def test_r_c_constructor
    c = Circle.new(r: 1.0, c: 0.0 + 0.0i)
    coeff = [1.0, 0.0, 0.0, -1.0]
    assert_equal c.coeff, coeff
    assert_equal c.r, 1.0
    assert_equal c.c, 0.0 + 0.0i
    c = Circle.new(r: nil, c: 1.0i)
    coeff = [0.0, -1.0+0.0i, -1.0+0.0i, 0.0]
    assert_equal c.coeff, coeff
    assert_nil c.r
    assert_equal c.c, 1.0i
  end
  
  def test_coeff_constructor
    coeff = [1.0, 0.0, 0.0, -1.0]
    c = Circle.new(coeff: coeff)
    assert_equal 1.0, c.r
    assert_equal 0.0 + 0.0i, c.c
    assert_equal coeff, c.coeff
    coeff = [0.0, -1.0+0.0i, -1.0+0.0i, 0.0]
    c = Circle.new(coeff: coeff)
    assert_nil c.r
    assert_equal 1.0i, c.c
    assert_equal coeff, c.coeff
  end

  def test_x_y_r_constructor
    c = Circle.new(x: 1.0, y: 2.0, r: 1.0)
    assert_equal 1.0, c.r
    assert_equal 1.0 + 2.0i, c.c
  end

  def test_is_geodesic
    c = Circle.new(x: Math::sqrt(2), y: 0.0, r: 1.0)
    assert c.is_geodesic?
    d = Circle.new(x: Math::sqrt(2), y: 0.0, r: 0.5)
    assert !d.is_geodesic?
    e = Circle.new(x: 0.5, y: 0.0, r: 1.0)
    assert !e.is_geodesic?
  end

  def test_in_disk
    p = 0.9 + 0.0i
    assert @c.in_disk?(p)
    p = -p
    assert !@c.in_disk?(p)
  end

  # So far only circles. Lines not implemented, yet.
  def test_triple_inclusion
    inner = Circle.new(x: 0.0, y: 0.0, r: 0.5)
    outer = Circle.new(x: 0.0, y: 0.5, r: 1.1)
    interdis = Circle.new(x: 0.0, y: 1.5, r: 0.5)
    p = 0.0 + 0.0i
    q = 0.0 + 1.0i
    r = 0.0 + 5.0i
    assert_equal 1, inner.triple_inclusion(p, outer)
    assert_equal 2, outer.triple_inclusion(p, inner)
    assert_equal 0, inner.triple_inclusion(p, interdis)
    assert_equal 0, outer.triple_inclusion(p, interdis)
    assert_equal 0, interdis.triple_inclusion(p, inner)
    assert_equal 2, inner.triple_inclusion(r, outer)
    assert_equal 1, outer.triple_inclusion(r, inner)
    assert_equal 0, inner.triple_inclusion(q, outer)
    assert_equal 0, outer.triple_inclusion(q, inner)
  end

  def test_translate_by
    # translation
    t = @prng.rand(-1000.0..1000.0)
    m = Matrix[[1.0, t], [0.0, 1.0]]
    c = @c.translate_by(m)
    assert_in_epsilon 1.0, c.r
    assert_in_epsilon @c.c + t, c.c
    # inversion
    m = Matrix[[0, 1], [1, 0]]
    c = @c.translate_by(m)
    4.times do |i|
      assert_in_epsilon @c.coeff[i], c.coeff[3-i]
    end
    # scaling
    m = Matrix[[t , 0.0], [0.0, 1.0 / t]]
    c = Circle.new(x: 0, y: 0, r: 1)
    c = c.translate_by(m)
    assert_in_epsilon 0.0, c.c
    assert_in_epsilon t**2, c.r
    # rotation
    m = Matrix[[0.0, 1.0], [-1.0, 0.0]]
    c = Circle.new(x: 0, y: 0, r: 1)
    d = c.translate_by(m)
    assert_in_epsilon c.r, d.r
    assert_in_epsilon c.c, d.c
    d = @c.translate_by(m)
    assert_in_epsilon @c.r, d.r
    assert_in_epsilon(-@c.c, d.c)
    x = exp(Math::PI / 4.0 * 1.0i)
    m = Matrix[[x, 0], [0, 1 / x]]
    d = c.translate_by(m)
    assert_in_epsilon c.r, d.r
    assert_in_epsilon(1i*c.c, d.c)
  end

  def test_get_bisector
    p = 0.5 + 0.0i
    q = -p
    c = Circle.get_bisector(p, q)
    assert_nil c.r
    assert_in_epsilon (p - q)*(-1i) / (p - q).abs, c.c
    p = 0.0 + 0.0i
    t = @prng.rand((0..(2*Math::PI)))
    q = Math::sqrt(0.5) * exp(t * 1.0i)
    c = Circle.get_bisector(p, q)
    assert_in_epsilon 1.0, c.r
    assert_in_epsilon Math::sqrt(2) * exp(t * 1.0i), c.c
  end
  
  def test_get_triple
    skip
  end

  def test_get_transformation
    r1 = @prng.rand(0.0..1000.0)
    r2 = @prng.rand(0.0..1000.0)
    t1 = @prng.rand(0.0..2*Math::PI)
    t2 = @prng.rand(0.0..2*Math::PI)
    c1 = Circle.new(c: Math::sqrt(1.0 + r1**2) * exp(t1), r: r1)
    c2 = Circle.new(c: Math::sqrt(1.0 + r2**2) * exp(t2), r: r2)
    c3 = Circle.new(c: exp(t2), r: nil)
    g = c1.get_transformation(c2)
    c = c1.translate_by(g)
    assert_in_epsilon c2.r, c.r
    assert_in_epsilon c2.c, c.c
    c = c2.translate_by(g.proj_inv)
    assert_in_epsilon c1.r, c.r
    assert_in_epsilon c1.c, c.c
    g = c1.get_transformation(c3)
    # c = c1.translate_by(g)
    # assert_nil c.r
    # assert_in_epsilon(-c3.c / c3.c.abs, c.c)
    c = c3.translate_by(g.proj_inv)
    assert_in_epsilon c1.r, c.r
    assert_in_epsilon c1.c, c.c
  end

  # TODO: Fix transformations
  def test_transformation_in_disk
    r1 = @prng.rand(0.0..10.0)
    r2 = @prng.rand(0.0..10.0)
    t1 = @prng.rand(0.0..2*Math::PI)
    t2 = @prng.rand(0.0..2*Math::PI)
    c1 = Circle.new(c: Math::sqrt(1.0 + r1**2) * exp(t1), r: r1)
    c2 = Circle.new(c: Math::sqrt(1.0 + r2**2) * exp(t2), r: r2)
    c3 = Circle.new(c: exp(t2), r: nil)
    g = c1.get_transformation(c2)
    100.times do |i|
      z = get_pr_in_disk
      w = g.moebius(z)
      assert w.abs < 1, "#{i} initial: #{z}, trafo: #{w}"
    end
    g = c1.get_transformation(c3)
    100.times do |i|
      z = get_pr_in_disk
      w = g.moebius(z)
      assert w.abs < 1, "#{i} initial: #{z}, trafo: #{w}"
    end
  end

  def get_pr_in_disk
    loop do
      x = @prng.rand(-1.0...1.0)
      y = @prng.rand(-1.0...1.0)
      z = x + y * 1.0i
      return z if z.abs < 1.0
    end
  end
  
  def test_plot_circle
    skip
  end
end
