require 'complex'

require 'daru'
include Daru

class Complex
  def plot_pt()
      pt = DataFrame.new({x: [real], y: [imag]})
      [pt, using: '2:3']
  end
end
