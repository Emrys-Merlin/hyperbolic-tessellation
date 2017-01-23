require_relative '../graph.rb'
require 'minitest/autorun'
require 'minitest/pride'

class TestGraph < MiniTest::Test
  def setup
    @p = 0.0 + 0.0i
    @gen = [Matrix[[(-2.220446049250313e-16+2.0i),
                    (5.551115123125783e-16-1.4142135623730954i)],
                   [(1.1102230246251565e-16-1.414213562373095i),
                    (-5.551115123125783e-17+1.9999999999999998i)]]]
    @graph = Graph.new(@p, @gen, 3)
  end

  def test_constructor
    graph = Graph.new(@p, @gen, 1)
    assert_equal @p, graph.root.p
    # This group is pretty symmetric, so we can check that the generators are
    # applied correctly by comparing the positions
    child1, child2 = graph.root.children
    assert_in_epsilon child1.p, -child2.p
    # deepen the graph to depth 2
    graph.deepen(2)
    # check that the matrix in every node really generates the orbit point
    graph.each do |node|
      assert_in_epsilon node.p, node.g.moebius(@p)
    end
  end

  def test_get_intersection
    @graph.get_intersection
    assert_equal 2, @graph.root.border.length
  end

  def test_get_tesselation
    @graph.get_tesselation
    @graph.each do |node|
      assert_equal 2, node.border.length
      @graph.root.border.length.times do |i|
        c = @graph.root.border[i]
        d = node.border[i]
        g = node.g
        c = c.translate_by(g)
        assert_in_epsilon d.r, c.r
        assert_in_epsilon d.c, c.c
      end
    end
  end
end
