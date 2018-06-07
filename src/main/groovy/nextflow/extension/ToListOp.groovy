/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.extension

import groovy.transform.CompileStatic
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.expression.DataflowExpression

/**
 * Implements {@link DataflowExtensions#toList(groovyx.gpars.dataflow.DataflowReadChannel)}  operator
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ToListOp {

    private DataflowReadChannel source

    private ordering

    /**
     * Implements both `toList` and `toSortedList`
     *
     * @param source The source channel
     * @param order Order the resulting channel, either boolean {@code true} or a comparing {@link Closure}
     */
    ToListOp( DataflowReadChannel source, order=null ) {
        this.source = source
        this.ordering = order
    }

    DataflowVariable apply() {
        assert source != null

        final target = new DataflowVariable()
        if( source instanceof DataflowExpression ) {
            final result = new ArrayList(1)
            Map<String,Closure> events = [:]
            events.onNext = { result.add(it) }
            events.onComplete = { target.bind(result) }
            DataflowHelper.subscribeImpl(source, events )
            return target
        }

        DataflowHelper.reduceImpl(source, target, []) { List list, item -> list << item }
        if( ordering ) {
            final sort = { List list -> ordering instanceof Closure ? list.sort((Closure) ordering) : list.sort() }
            return (DataflowVariable)target.then(sort)
        }
        else {
            return target
        }
    }

    static DataflowVariable apply( DataflowReadChannel source ) {
        new ToListOp(source).apply()
    }

}
