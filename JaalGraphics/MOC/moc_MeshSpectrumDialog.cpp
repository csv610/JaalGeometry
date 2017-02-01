/****************************************************************************
** Meta object code from reading C++ file 'MeshSpectrumDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshSpectrumDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshSpectrumDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshSpectrumDialog_t {
    QByteArrayData data[4];
    char stringdata0[43];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshSpectrumDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshSpectrumDialog_t qt_meta_stringdata_JMeshSpectrumDialog = {
    {
        QT_MOC_LITERAL(0, 0, 19), // "JMeshSpectrumDialog"
        QT_MOC_LITERAL(1, 20, 10), // "genVectors"
        QT_MOC_LITERAL(2, 31, 0), // ""
        QT_MOC_LITERAL(3, 32, 10) // "viewVector"

    },
    "JMeshSpectrumDialog\0genVectors\0\0"
    "viewVector"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshSpectrumDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    2,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   24,    2, 0x08 /* Private */,
    3,    0,   25,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JMeshSpectrumDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshSpectrumDialog *_t = static_cast<JMeshSpectrumDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->genVectors();
            break;
        case 1:
            _t->viewVector();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshSpectrumDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JMeshSpectrumDialog.data,
        qt_meta_data_JMeshSpectrumDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JMeshSpectrumDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshSpectrumDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshSpectrumDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshSpectrumDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshSpectrumDialog"))
        return static_cast< Ui::MeshSpectrumDialog*>(const_cast< JMeshSpectrumDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshSpectrumDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 2)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 2;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
