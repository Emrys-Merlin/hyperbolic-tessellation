# frozen_string_literal: true
require_relative 'circle.rb'
require_relative 'graph.rb'

require 'gnuplotrb'
include GnuplotRB

# Plot class
class FDPlot
  attr_reader :plot, :rplot, :graph

  def initialize(p, gen, depth = 3, points = true)
    @graph = Graph.new(p, gen, depth)
    @graph.get_tessellation

    @plot = plot_boundary

    @graph.each do |node|
      node.border.each do |circ|
        @plot = @plot << circ.plot_circle
      end
    end

    plot_points if points

    @rplot = plot_fd
  end

  def plot_fd(center = @graph.root)
    plot = plot_boundary
    center.border.each do |circ|
      plot = plot << circ.plot_circle
    end
    plot
  end

  private

  def plot_boundary
    x = []
    y = []
    (0..1000).each do |j|
      z = exp(2i * Math::PI / 1000 * j)
      x << z.real
      y << z.imag
    end
    @boundary = [x, y]

    @df = [@boundary, with: 'lines', lt: "rgb 'black'"]
    Plot.new(@df,
             xrange: -1..1,
             yrange: -1..1,
             border: '0',
             ytics: nil,
             xtics: nil,
             key: nil)
  end

  def plot_points
    @graph.each do |node|
      @plot = @plot << plot_pt(node.p)
    end
  end

  def plot_pt(p)
    df = [[p.real], [p.imag]]
    [df]
  end
end
