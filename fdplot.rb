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
        # TODO: Fix large radii (best if they could be prevented)
        next if circ.r > 1e6
        @plot = @plot << circ.plot_circle(node.angles[circ])
      end
    end

    plot_points if points

    @rplot = plot_fd
  end

  def plot_fd(center = @graph.root)
    plot = plot_boundary
    center.border.each do |circ|
      plot = plot << circ.plot_circle(center.angles[circ])
    end
    plot
  end

  def plot_hemisphere
    plot = plot_boundary_hs
    @graph.each do |node|
      node.border.each do |circ|
        plot = plot << circ.plot_circle_hs(node.angles[circ])
      end
    end
    plot
  end

  def plot_hemisphere_rhino(filename, step = 1000, scale = 10, eps = 0.001)
    File.open(filename, 'w') do |f|
      f.puts('_InterpCrv')
      (0..step).each do |j|
        w = scale * exp(2i * Math::PI / step * j)
        f.puts("#{w.real},#{w.imag},0")
      end
      f.puts('_Enter')
      @graph.each do |node|
        circles = []
        same = false
        node.border.each do |circ|
          circles.each do |c|
            if (c.r - circ.r).abs < eps && (c.c - circ.c).abs < eps
              same = true
              break
            end
          end
          next if same
          circles << circ
          f.puts('_InterpCrv')
          circ.generalized_circle_points_hs(node.angles[circ],step).transpose.each do |pt|
            f.puts("#{scale * pt[0]},#{scale * pt[1]},#{scale * pt[2]}")
          end
          f.puts('_Enter')
        end
      end
    end
  end
  
  private

  def plot_boundary
    x = []
    y = []
    (0...1000).each do |j|
      z = exp(2i * Math::PI / 1000 * j)
      x << z.real
      y << z.imag
    end

    x[-1] = 1
    y[-1] = 0
    
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

  def plot_boundary_hs
    x = []
    y = []
    z = []
    (0...1000).each do |j|
      w = exp(2i * Math::PI / 1000 * j)
      x << w.real
      y << w.imag
      z << 0
    end
    @boundary = [x, y, z]

    @df = [@boundary, with: 'lines', lt: "rgb 'black'"]
    Splot.new(@df,
             xrange: -1..1,
             yrange: -1..1,
             zrange: -1..1,
             style: 'data lines',
             hidden3d: true,
             border: '0',
             ytics: nil,
             xtics: nil,
             ztics: nil,
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
