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

package nextflow.mail
import java.nio.file.Path

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.util.CheckHelper

/**
 * Helper class modeling mail parameters
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class Mail {

    private static final Map ATTACH_HEADERS = [
            contentId: String,
            disposition:String,
            fileName: String,
            description: String
    ]

    String from

    String to

    String cc

    String bcc

    String subject

    String charset

    String body

    String text

    List<Attachment> attachments

    String type

    /**
     * Creates a {@link Mail} object given a {@link Mail} object
     *
     * @param params
     * @return A mail object representing the message to send
     */
    static Mail of(Map params) {
        def result = new Mail()

        if( params.from )
            result.from(params.from.toString())

        if( params.to )
            result.to(params.to.toString())

        if( params.cc )
            result.cc(params.cc.toString())

        if( params.bcc )
            result.bcc(params.bcc.toString())

        if( params.subject )
            result.subject(params.subject.toString())

        if( params.charset )
            result.charset(params.charset.toString())

        if( params.type )
            result.type(params.type.toString())

        if( params.body )
            result.body(params.body)

        if( params.text )
            result.text(params.text)

        if( params.attach )
            result.attach(params.attach)

        return result
    }

    /**
     * Mail sender (addresses must follow RFC822 syntax)
     * Multiple addresses can be separated by a comma
     */
    void from( String address ) {
        this.from = address
    }

    /**
     * Mail TO recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    void to( String address ) {
        this.to = address
    }

    /**
     * Mail CC recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    void cc( String address ) {
        this.cc = address
    }

    /**
     * Mail BCC recipients (addresses must follow RFC822 syntax).
     * Multiple addresses can be separated by a comma.
     */
    void bcc( String address ) {
        this.bcc = address
    }

    /**
     * @param subject The email subject
     */
    void subject( String subject ) {
        this.subject = subject
    }

    private String stringify( value ) {
        if( value instanceof File )
            return value.text
        if( value instanceof Path )
            return value.text
        if( value instanceof CharSequence )
            return value.toString()
        if( value != null )
            throw new IllegalArgumentException("Not a valid mail body argument [${value.getClass().getName()}]: $value")
        return null
    }

    /**
     * @param str The email content
     */
    void body( value ) {
        this.body = stringify(value)
    }

    /**
     * Plain text mail content
     * @param text The mail text content
     */
    void text( value ) {
        this.text = stringify(value)
    }

    /**
     * Mail content mime-type
     */
    void type( String mime )  {
        this.type = mime
    }

    /**
     * @param charset The mail content charset
     */
    void charset( String charset ) {
        this.charset = charset
    }

    /**
     * Add an email attachment
     *
     * @param item A attachment file path either as {@link File}, {@code Path} or {@link String} path
     */
    void attach( item ) {
        if( this.attachments == null )
            this.attachments = []

        if( item instanceof Attachment ) {
            this.attachments << (Attachment)item
        }
        else if( item instanceof Collection ) {
            this.attachments.addAll( item.collect{ new Attachment(it) } )
        }
        else if( item instanceof Object[] ) {
            this.attachments.addAll( item.collect{ new Attachment(it) } )
        }
        else if( item ) {
            this.attachments << new Attachment(item)
        }
    }

    /**
     * Add an email attachment headers
     *
     * @param headers
     *      Attachment optional content directives. The following parameters are accepted:
     *      - contentId:  Set the "Content-ID" header field of this body part
     *      - fileName:  Set the filename associated with this body part, if possible
     *      - description: Set the "Content-Description" header field for this body part
     *      - disposition: Set the "Content-Disposition" header field of this body part
     *
     * @param item
     */
    void attach( Map headers, item ) {
        CheckHelper.checkParams('attach', headers, ATTACH_HEADERS)

        if( this.attachments == null )
            this.attachments = []

        this.attachments << new Attachment(headers, item)
    }

}
