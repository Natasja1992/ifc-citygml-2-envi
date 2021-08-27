from conversionTool import ConversionTool
from emlWriter import EmlWriter
import argparse


def main():
    ###################################################################################################################
    # User input
    ###################################################################################################################
    parser = argparse.ArgumentParser()

    parser.add_argument("--ifc", help="Give path to input IFC file", required=True)
    parser.add_argument("--citygml", help="Give path to input CityGMl file", required=True)
    parser.add_argument("--output", help="Give path to output file", required=True)
    parser.add_argument("--borderwidth", help="Set borderwidth around site model in meters (default: 20.0)", default=20.0)
    parser.add_argument("--resolution", help="Set resolution of the output model grid (default: 1.0)", default=1.0)
    parser.add_argument("--bordergrid", help="Set number of empty grid cells around the model (default: 5)", default=5)

    args = parser.parse_args()

    variables = {
        "citygml_borderwidth": float(args.borderwidth),
        "resolution": (float(args.resolution), float(args.resolution), float(args.resolution)),
        "border_grid": int(args.bordergrid)
    }

    ###################################################################################################################

    conversion_tool = ConversionTool(args.ifc, args.citygml, variables)
    eml_writer = EmlWriter(conversion_tool)
    eml_writer.write(args.output)


if __name__ == '__main__':
    main()
