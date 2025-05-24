import React, {useContext, useEffect, useRef, useState} from 'react'
import {useParams} from "react-router-dom";
import {Spinner} from '@radix-ui/themes'
import {AppContext} from "../App";

export default function ToolPage() {
    const {toolURL} = useParams();

    const toolsList = useContext(AppContext);
    console.log(toolsList)

    const [currentTool, setCurrentTool] = useState(null);
    const [inputText, setInputText] = useState(null);
    const [outputText, setOutputText] = useState('');
    const [isLoading, setIsLoading] = useState(false);

    const textInputRef = useRef(null);
    const fileInputRef = useRef(null);

    //initializes currentTool
    useEffect(() => {
        toolsList.forEach((tool) => {
            if (tool.URL === toolURL) {
                setCurrentTool(tool);
            }
        });
    }, []);

    //prevents page from rendering until after the currentTool state is set
    if (!currentTool) return;



    function handleTextChange(e) {
        setInputText(e.target.value);
        //logic for only allowing one input type at a time
        if (e.target.value === '') {
            fileInputRef.current.removeAttribute('disabled');
        }
        else {
            fileInputRef.current.setAttribute('disabled', 'true');
        }
    }



    function handleFileChange(event) {
        const reader = new FileReader();
        try {
            reader.onload = (e) => {
                setInputText(e.target.result);
            }
            reader.readAsText(event.target.files[0]);
        }
        catch (error) {
            setOutputText(error);
        }
        //logic for only allowing one input type at a time
        if (event.target.value === '') {
            textInputRef.current.removeAttribute('disabled');
        }
        else {
            textInputRef.current.setAttribute('disabled', 'true');
        }
    }



    function handleRemoveFile() {
        setInputText('');
        fileInputRef.current.value = '';
        textInputRef.current.removeAttribute('disabled');
    }



    function handleSubmit(e) {
        e.preventDefault();
        setIsLoading(true);

        fetch(`http://localhost:5000/tool-page/submit/${currentTool.URL}`, {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({'inputText': inputText})
        })
            .then((response) => response.json())
            .then((data) => {
                setOutputText(data.result);
                setIsLoading(false);
            })
            .catch((error) => {
                console.log(error.message)
                setOutputText('error fetching algorithm from server');
                setIsLoading(false);
            })
    }



    return (
        <section className={'container tool-page'}>
            <h1>{currentTool.toolTitle}</h1>
            <p>{currentTool.longDescription}</p>
            <div className={'container tool-form-container'}>
                <p>Input format:</p>
                <p className={'input-format'}>{currentTool.inputFormat}</p>
                <p>Output format:</p>
                <p className={'output-format'}>{currentTool.outputFormat}</p>
                <p>Example:</p>
                <code className={'input-example'}>{currentTool.inputExample}</code>
                <p className={'input-instructions'}>Input either in text box or as .txt file</p>
                <hr/>
                <form className={'container tool-form'} onSubmit={handleSubmit}>
                    <textarea className={'algorithm-text-input'} ref={textInputRef} onChange={handleTextChange}/>
                    <input className={'algorithms-file-input'} ref={fileInputRef} type={'file'} onChange={handleFileChange}/>
                    <button className={'remove-file-button'} type={'button'} onClick={handleRemoveFile}>Remove file</button>
                    <button className={'blue-button'} type={'submit'}>Calculate</button>
                </form>
            </div>
            {isLoading ? <Spinner /> : ''}
            <p>Result:</p>
            <p className={'output-text'}>{outputText}</p>
        </section>
    );
}