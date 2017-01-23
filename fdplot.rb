require_relative 'circle.rb'
require_relative 'graph.rb'

require_relative 'complex.rb'

require 'gnuplotrb'
include GnuplotRB

class FDPlot < Plot
  def initialize(p, gen, depth = 3)
    @graph = Graph.new(p, gen, depth)
    @graph.get_tesselation

    @graph.each do |node|
      
    end
    
    super(plot_boundary,
          xrange: -1..1,
           yrange: -1..1,
           border: '0',
           ytics: nil,
           xtics: nil)
  end

  private
  def plot_boundary
    @boundary = Circle.new(x: 0, y: 0, r: 1)
    @boundary.plot_circle('black')
  end
end
